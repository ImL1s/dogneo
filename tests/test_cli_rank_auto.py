"""Tests for `dogneo rank` auto-detection features.

TDD RED: Tests for automatic proteome/alleles detection in `rank` command.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from dogneo.cli import cli


class TestRankAutoDetection:
    """Tests for rank command's auto-detection behavior."""

    def test_rank_auto_proteome_from_cache(self, tmp_path):
        """rank without --protein-db auto-detects cached proteome."""
        # Set up fake cache with proteome
        cache = tmp_path / "cache"
        cache.mkdir()
        fa = cache / "CanFam3.1.pep.all.fa"
        # Write a minimal but valid proteome with Ensembl-style headers
        fa.write_text(
            ">ENSCAFP00000000001 gene_symbol:TP53\n"
            "MKTAYIAKQRQISFVKSHFSRQ\n"
        )

        demo_dir = Path(__file__).parent.parent / "dogneo" / "data" / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', cache)):
            result = runner.invoke(cli, [
                "rank",
                "--vcf", str(vcf_path),
                "--output-dir", str(tmp_path / "out"),
            ])

        # Should succeed and mention auto-loading proteome
        assert result.exit_code == 0
        assert "auto" in result.output.lower() or "cached" in result.output.lower() or "Loading protein" in result.output

    def test_rank_auto_alleles_from_bundled(self, tmp_path):
        """rank without --alleles uses bundled DLA alleles."""
        demo_dir = Path(__file__).parent.parent / "dogneo" / "data" / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"
        proteome = Path(__file__).parent.parent / "data" / "reference" / "CanFam3.1.pep.all.fa"

        if not vcf_path.exists() or not proteome.exists():
            pytest.skip("Demo data not available")

        runner = CliRunner()
        result = runner.invoke(cli, [
            "rank",
            "--vcf", str(vcf_path),
            "--protein-db", str(proteome),
            "--output-dir", str(tmp_path / "out"),
        ])

        # Should succeed even without --alleles
        assert result.exit_code == 0
        assert "DLA" in result.output or "allele" in result.output.lower()

    def test_rank_no_proteome_no_cache_warns(self, tmp_path):
        """rank without proteome nor cache shows helpful message."""
        demo_dir = Path(__file__).parent.parent / "dogneo" / "data" / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        runner = CliRunner()
        with patch("dogneo.data.manager.ReferenceDataManager.__init__",
                   lambda self, cache_dir=None: setattr(self, 'cache_dir', tmp_path / "empty")):
            result = runner.invoke(cli, [
                "rank",
                "--vcf", str(vcf_path),
                "--output-dir", str(tmp_path / "out"),
            ])

        # Should still succeed (no peptides) but mention setup
        assert "setup" in result.output.lower() or "protein" in result.output.lower()

    def test_rank_explicit_protein_db_overrides_cache(self, tmp_path):
        """Explicit --protein-db still works and takes precedence."""
        demo_dir = Path(__file__).parent.parent / "dogneo" / "data" / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"
        proteome = Path(__file__).parent.parent / "data" / "reference" / "CanFam3.1.pep.all.fa"

        if not vcf_path.exists() or not proteome.exists():
            pytest.skip("Demo data not available")

        runner = CliRunner()
        result = runner.invoke(cli, [
            "rank",
            "--vcf", str(vcf_path),
            "--protein-db", str(proteome),
            "--output-dir", str(tmp_path / "out"),
        ])

        assert result.exit_code == 0
        assert "candidates" in result.output.lower() or "✅" in result.output
