"""Tests for dogneo.llm.cli_wrapper — CLI command building and invocation."""
from __future__ import annotations

import subprocess

from dogneo.llm.cli_wrapper import (
    _build_claude_cmd,
    _build_codex_cmd,
    _build_env,
    _build_gemini_cmd,
    call_ai_cli,
    call_ai_with_fallback,
    check_cli_availability,
)

# ---------------------------------------------------------------------------
# Command building
# ---------------------------------------------------------------------------

class TestBuildGeminiCmd:
    """Tests for Gemini CLI command construction."""

    def test_default_model(self):
        cmd = _build_gemini_cmd("hello")
        assert cmd[0] == "gemini"
        assert "-p" in cmd
        assert "hello" in cmd
        assert "gemini-2.5-flash" in cmd

    def test_custom_model(self):
        cmd = _build_gemini_cmd("hello", model="gemini-2.5-pro")
        assert "gemini-2.5-pro" in cmd

    def test_json_output(self):
        cmd = _build_gemini_cmd("hello", output_format="json")
        assert "json" in cmd


class TestBuildClaudeCmd:
    """Tests for Claude Code CLI command construction."""

    def test_default(self):
        cmd = _build_claude_cmd("test prompt")
        assert cmd[0] == "claude"
        assert "-p" in cmd
        assert "test prompt" in cmd
        assert "--max-turns" in cmd
        assert "1" in cmd

    def test_output_format(self):
        cmd = _build_claude_cmd("test", output_format="json")
        assert "--output-format" in cmd
        idx = cmd.index("--output-format")
        assert cmd[idx + 1] == "json"


class TestBuildCodexCmd:
    """Tests for Codex CLI command construction."""

    def test_default(self):
        cmd = _build_codex_cmd("test prompt")
        assert cmd[0] == "codex"
        assert "exec" in cmd
        assert "--skip-git-repo-check" in cmd
        assert "-s" in cmd

    def test_custom_model(self):
        cmd = _build_codex_cmd("hello", model="o4-mini")
        assert "o4-mini" in cmd


# ---------------------------------------------------------------------------
# Environment building
# ---------------------------------------------------------------------------

class TestBuildEnv:
    """Tests for PATH environment construction."""

    def test_env_has_path(self):
        env = _build_env()
        assert "PATH" in env

    def test_env_includes_homebrew(self):
        env = _build_env()
        # At least one of these should be in PATH
        assert any(p in env["PATH"] for p in ["/opt/homebrew/bin", "/usr/local/bin"])


# ---------------------------------------------------------------------------
# call_ai_cli with monkeypatch
# ---------------------------------------------------------------------------

class TestCallAiCli:
    """Tests for call_ai_cli subprocess invocation."""

    def test_success(self, monkeypatch):
        """Successful CLI call returns output."""
        def mock_run(cmd, **kwargs):
            return subprocess.CompletedProcess(
                args=cmd, returncode=0,
                stdout="Generated answer", stderr="",
            )
        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_cli("test prompt", tool="gemini")
        assert result.success is True
        assert result.output == "Generated answer"
        assert result.tool == "gemini"

    def test_nonzero_exit(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            return subprocess.CompletedProcess(
                args=cmd, returncode=1,
                stdout="", stderr="Auth error",
            )
        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_cli("test", tool="claude")
        assert result.success is False
        assert "Auth error" in result.error

    def test_timeout(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            raise subprocess.TimeoutExpired(cmd=cmd, timeout=10)
        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_cli("test", tool="gemini", timeout=10)
        assert result.success is False
        assert "timed out" in result.error

    def test_file_not_found(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            raise FileNotFoundError("gemini not found")
        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_cli("test", tool="gemini")
        assert result.success is False
        assert "not found" in result.error

    def test_empty_output(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            return subprocess.CompletedProcess(
                args=cmd, returncode=0, stdout="", stderr="",
            )
        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_cli("test", tool="codex")
        assert result.success is False
        assert "Empty output" in result.error

    def test_unknown_tool(self):
        result = call_ai_cli("test", tool="unknown_tool_xyz")
        assert result.success is False
        assert "Unknown CLI tool" in result.error


# ---------------------------------------------------------------------------
# call_ai_with_fallback
# ---------------------------------------------------------------------------

class TestCallAiWithFallback:
    """Tests for the CLI fallback chain."""

    def test_first_tool_succeeds(self, monkeypatch):
        call_count = 0

        def mock_run(cmd, **kwargs):
            nonlocal call_count
            call_count += 1
            return subprocess.CompletedProcess(
                args=cmd, returncode=0,
                stdout="gemini answer", stderr="",
            )

        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_with_fallback("prompt")
        assert result.success is True
        assert result.output == "gemini answer"
        assert call_count == 1  # Only tried first tool

    def test_fallback_to_second(self, monkeypatch):
        attempts = []

        def mock_run(cmd, **kwargs):
            tool_name = cmd[0]
            attempts.append(tool_name)
            if tool_name == "gemini":
                raise FileNotFoundError("not found")
            return subprocess.CompletedProcess(
                args=cmd, returncode=0,
                stdout="claude answer", stderr="",
            )

        monkeypatch.setattr(subprocess, "run", mock_run)
        result = call_ai_with_fallback("prompt")
        assert result.success is True
        assert "claude answer" in result.output
        assert "gemini" in attempts
        assert "claude" in attempts


# ---------------------------------------------------------------------------
# check_cli_availability
# ---------------------------------------------------------------------------

class TestCheckCliAvailability:
    """Tests for CLI tool availability checking."""

    def test_all_available(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            return subprocess.CompletedProcess(
                args=cmd, returncode=0, stdout="v1.0.0", stderr="",
            )
        monkeypatch.setattr(subprocess, "run", mock_run)
        avail = check_cli_availability()
        assert avail["gemini"] is True
        assert avail["claude"] is True
        assert avail["codex"] is True

    def test_none_available(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            raise FileNotFoundError("not found")
        monkeypatch.setattr(subprocess, "run", mock_run)
        avail = check_cli_availability()
        assert avail["gemini"] is False
        assert avail["claude"] is False
        assert avail["codex"] is False

    def test_partial_availability(self, monkeypatch):
        def mock_run(cmd, **kwargs):
            if cmd[0] == "gemini":
                return subprocess.CompletedProcess(
                    args=cmd, returncode=0, stdout="v1.0", stderr="",
                )
            raise FileNotFoundError("not found")
        monkeypatch.setattr(subprocess, "run", mock_run)
        avail = check_cli_availability()
        assert avail["gemini"] is True
        assert avail["claude"] is False
        assert avail["codex"] is False
