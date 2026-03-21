"""Subprocess wrapper for AI CLI tools (Gemini CLI, Claude Code, Codex).

Provides free, no-API-key LLM access via official CLI tools running in
headless/non-interactive mode. Handles PATH setup for cron/SSH environments,
timeout management, and fallback chains.

Reference: cli-as-ai-backend pattern — uses official CLI flags only (TOS-compliant).
"""

from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class CLIResult:
    """Result from an AI CLI invocation.

    Attributes:
        output: Generated text output.
        tool: CLI tool name used.
        model: Model name used.
        success: Whether the call succeeded.
        error: Error message if failed.
    """

    output: str = ""
    tool: str = ""
    model: str = ""
    success: bool = True
    error: str = ""


def _build_env() -> dict[str, str]:
    """Build environment with PATH that includes NVM node binaries.

    This is critical for cron/SSH environments where interactive
    shell paths are not available.
    """
    env = os.environ.copy()
    home = os.environ.get("HOME", "")

    # Common NVM node binary paths
    nvm_paths = [
        f"{home}/.nvm/versions/node/v22.20.0/bin",
        f"{home}/.nvm/versions/node/v20.18.0/bin",
        f"{home}/.nvm/versions/node/v18.20.0/bin",
    ]

    # Homebrew paths
    brew_paths = [
        "/opt/homebrew/bin",
        "/usr/local/bin",
    ]

    existing_path = env.get("PATH", "/usr/bin:/bin")
    extra_paths = [p for p in nvm_paths + brew_paths if p not in existing_path]

    if extra_paths:
        env["PATH"] = ":".join(extra_paths) + ":" + existing_path

    return env


def _build_gemini_cmd(
    prompt: str,
    model: str = "gemini-2.5-flash",
    output_format: str = "text",
) -> list[str]:
    """Build Gemini CLI command.

    Flags:
        -m: model selection
        -p: non-interactive prompt (headless mode)
        -o: output format (text/json/stream-json)
    """
    return ["gemini", "-m", model, "-p", prompt, "-o", output_format]


def _build_claude_cmd(
    prompt: str,
    model: str = "claude-sonnet-4-6",
    output_format: str = "text",
) -> list[str]:
    """Build Claude Code CLI command.

    Flags:
        -p: non-interactive prompt (headless mode)
        --model: model selection
        --output-format: text/json/stream-json
        --max-turns 1: limit to single turn for non-agentic use
    """
    return [
        "claude", "-p", prompt,
        "--model", model,
        "--output-format", output_format,
        "--max-turns", "1",
    ]


def _build_codex_cmd(
    prompt: str,
    model: str = "gpt-5.4",
) -> list[str]:
    """Build Codex CLI command.

    Flags:
        exec: non-interactive mode
        -m: model selection
        --skip-git-repo-check: run outside git repo
        -s read-only: safe sandbox mode
    """
    return [
        "codex", "exec", prompt,
        "-m", model,
        "--skip-git-repo-check",
        "-s", "read-only",
    ]


def call_ai_cli(
    prompt: str,
    tool: str = "gemini",
    model: str | None = None,
    timeout: int = 180,
    output_format: str = "text",
) -> CLIResult:
    """Call an AI CLI tool via subprocess.

    This uses official CLI tools in their documented headless mode,
    which is fully TOS-compliant.

    Args:
        prompt: The prompt text to send.
        tool: CLI tool name ("gemini", "claude", "codex").
        model: Model to use. Defaults vary by tool.
        timeout: Maximum seconds to wait.
        output_format: Output format ("text" or "json").

    Returns:
        CLIResult with output or error.
    """
    # Default models per tool
    default_models = {
        "gemini": "gemini-2.5-flash",
        "claude": "claude-sonnet-4-6",
        "codex": "gpt-5.4",
    }

    if model is None:
        model = default_models.get(tool, "")

    # Build command
    if tool == "gemini":
        cmd = _build_gemini_cmd(prompt, model, output_format)
    elif tool == "claude":
        cmd = _build_claude_cmd(prompt, model, output_format)
    elif tool == "codex":
        cmd = _build_codex_cmd(prompt, model)
    else:
        return CLIResult(
            success=False,
            error=f"Unknown CLI tool: {tool}",
            tool=tool,
            model=model,
        )

    env = _build_env()

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )

        # stderr may contain status messages (e.g., "Loaded cached credentials")
        # Only treat as error if return code is non-zero
        if result.returncode != 0:
            error_msg = result.stderr.strip() or f"Exit code: {result.returncode}"
            logger.warning("%s CLI error: %s", tool, error_msg)
            return CLIResult(
                success=False,
                error=error_msg,
                tool=tool,
                model=model,
            )

        output = result.stdout.strip()
        if not output:
            return CLIResult(
                success=False,
                error="Empty output from CLI",
                tool=tool,
                model=model,
            )

        return CLIResult(
            output=output,
            tool=tool,
            model=model,
            success=True,
        )

    except FileNotFoundError:
        msg = f"{tool} CLI not found. Install: "
        install_cmds = {
            "gemini": "npm install -g @google/gemini-cli",
            "claude": "npm install -g @anthropic-ai/claude-code",
            "codex": "npm install -g @openai/codex",
        }
        msg += install_cmds.get(tool, f"see {tool} documentation")
        logger.warning(msg)
        return CLIResult(success=False, error=msg, tool=tool, model=model)

    except subprocess.TimeoutExpired:
        msg = f"{tool} CLI timed out after {timeout}s"
        logger.warning(msg)
        return CLIResult(success=False, error=msg, tool=tool, model=model)


def call_ai_with_fallback(
    prompt: str,
    tools_chain: list[tuple[str, str]] | None = None,
    timeout: int = 180,
    output_format: str = "text",
) -> CLIResult:
    """Try multiple AI CLIs in order, returning first successful result.

    Args:
        prompt: The prompt text.
        tools_chain: List of (tool_name, model_name) tuples to try in order.
            Defaults to Gemini → Claude → Codex.
        timeout: Timeout per attempt.
        output_format: Output format.

    Returns:
        CLIResult from the first successful tool, or last failure.
    """
    if tools_chain is None:
        tools_chain = [
            ("gemini", "gemini-2.5-flash"),
            ("claude", "claude-sonnet-4-6"),
            ("codex", "gpt-5.4"),
        ]

    last_result = CLIResult(success=False, error="No tools configured")

    for tool_name, model_name in tools_chain:
        logger.debug("Trying %s (%s)...", tool_name, model_name)
        result = call_ai_cli(
            prompt=prompt,
            tool=tool_name,
            model=model_name,
            timeout=timeout,
            output_format=output_format,
        )

        if result.success:
            logger.info("Success with %s (%s)", tool_name, model_name)
            return result

        logger.debug("%s failed: %s", tool_name, result.error)
        last_result = result

    logger.error("All CLI tools failed. Last error: %s", last_result.error)
    return last_result


def check_cli_availability() -> dict[str, bool]:
    """Check which AI CLI tools are available on this system.

    Returns:
        Dict mapping tool name to availability.
    """
    env = _build_env()
    available: dict[str, bool] = {}

    checks = {
        "gemini": ["gemini", "--version"],
        "claude": ["claude", "--version"],
        "codex": ["codex", "--version"],
    }

    for tool, cmd in checks.items():
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=10, env=env,
            )
            available[tool] = result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            available[tool] = False

    return available
