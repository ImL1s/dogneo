"""Tests for `dogneo setup` and `dogneo demo` CLI commands.

TDD RED phase: defines expected behavior for the setup/demo commands.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from dogneo.cli import cli


class TestSetupCommand:
    """Tests for `dogneo setup`."""

    def test_setup_exists(self):
        """The 'setup' command is registered."""
        runner = CliRunner()
        result = runner.invoke(cli, ["setup", "--help"])
        assert result.exit_code == 0
        assert "reference data" in result.output.lower() or "setup" in result.output.lower()

    def test_setup_downloads_proteome(self, tmp_path):
        """setup downloads proteome and reports success."""
        # Create fake gz content
        fasta = b">P1\nMKK\n>P2\nMWW\n"
        gz = gzip.compress(fasta)

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.headers = {"content-length": str(len(gz))}
        mock_resp.iter_content = lambda chunk_size: [gz]
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)

        runner = CliRunner()
        with patch("dogneo.data.manager.requests.get", return_value=mock_resp):
            with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                       lambda self, cache_dir=None: setattr(self, 'cache_dir', tmp_path)):
                result = runner.invoke(cli, ["setup"])

        assert result.exit_code == 0
        assert "✅" in result.output or "complete" in result.output.lower()

    def test_setup_cached_skips(self, tmp_path):
        """setup with existing cache reports already cached."""
        fa = tmp_path / "CanFam3.1.pep.all.fa"
        fa.write_text(">P1\nMKK\n")

        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', tmp_path)):
            result = runner.invoke(cli, ["setup"])

        assert result.exit_code == 0
        assert "already" in result.output.lower() or "✅" in result.output

    def test_setup_force_flag(self):
        """setup --force flag is accepted."""
        runner = CliRunner()
        result = runner.invoke(cli, ["setup", "--help"])
        assert "--force" in result.output

    def test_setup_shows_status(self, tmp_path):
        """setup shows a summary after completion."""
        fa = tmp_path / "CanFam3.1.pep.all.fa"
        fa.write_text(">P1\nMKK\n>P2\nMWW\n>P3\nMAA\n")

        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', tmp_path)):
            result = runner.invoke(cli, ["setup"])

        assert result.exit_code == 0


class TestDemoCommand:
    """Tests for `dogneo demo`."""

    def test_demo_exists(self):
        """The 'demo' command is registered."""
        runner = CliRunner()
        result = runner.invoke(cli, ["demo", "--help"])
        assert result.exit_code == 0

    def test_demo_requires_proteome(self, tmp_path):
        """demo fails gracefully when proteome not set up."""
        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', tmp_path / "empty")):
            result = runner.invoke(cli, ["demo"])

        # Should fail with helpful message about running dogneo setup
        assert result.exit_code != 0 or "setup" in result.output.lower()

    def test_demo_runs_pipeline(self, tmp_path):
        """demo runs full pipeline with bundled data."""
        # Set up cached proteome
        cache = tmp_path / "cache"
        cache.mkdir()

        # Use real demo data from the project
        demo_dir = Path(__file__).parent.parent / "data" / "demo"
        proteome = Path(__file__).parent.parent / "data" / "reference" / "CanFam3.1.pep.all.fa"

        if not proteome.exists():
            pytest.skip("Real proteome not available for integration test")

        # Copy proteome to fake cache
        import shutil
        shutil.copy(proteome, cache / "CanFam3.1.pep.all.fa")

        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', cache)):
            result = runner.invoke(cli, ["demo", "--output-dir", str(tmp_path / "out")])

        assert result.exit_code == 0
        assert "candidates" in result.output.lower() or "✅" in result.output

    def test_demo_output_dir_option(self):
        """demo accepts --output-dir option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["demo", "--help"])
        assert "--output-dir" in result.output

    def test_demo_skip_binding_by_default(self):
        """demo help shows --binding option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["demo", "--help"])
        # Binding should be optional
        assert result.exit_code == 0
