"""Tests for demo cache behavior.

TDD Phase 1.2: demo should reuse cached proteome, not re-download.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest


class TestDemoCacheBehavior:
    """Ensure demo reuses cached proteome and does not re-download."""

    def test_demo_skips_download_when_cached(self, tmp_path):
        """When proteome is already cached, demo must not trigger download."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager()
        cached = mgr.get_proteome_path()
        if not cached:
            pytest.skip("Proteome not cached — run dogneo setup first")

        from click.testing import CliRunner
        from dogneo.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["demo", "--output-dir", str(tmp_path / "demo_out")])

        # Output should NOT contain download messages
        assert "Downloading" not in result.output
        assert "Decompressing" not in result.output

    def test_demo_shows_using_cached_message(self, tmp_path):
        """Demo should indicate it's using cached proteome."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager()
        cached = mgr.get_proteome_path()
        if not cached:
            pytest.skip("Proteome not cached — run dogneo setup first")

        from click.testing import CliRunner
        from dogneo.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["demo", "--output-dir", str(tmp_path / "demo_out")])

        # Should show cached/auto-detected message
        assert "Using" in result.output or "cached" in result.output.lower() or "Auto-detected" in result.output
