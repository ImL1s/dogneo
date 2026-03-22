"""Tests for the rank pipeline service layer.

TDD Phase 1.1: Extract rank logic from cli.py into a reusable service.
"""
from __future__ import annotations

from pathlib import Path

import pytest


class TestRankInput:
    """Test RankInput dataclass validation."""

    def test_rank_input_has_required_fields(self):
        from dogneo.app.rank_pipeline import RankInput

        inp = RankInput(vcf_path=Path("test.vcf"))
        assert inp.vcf_path == Path("test.vcf")
        assert inp.sample_id == "SAMPLE"
        assert inp.alleles == []
        assert inp.mhci_lengths == [8, 9, 10, 11]
        assert inp.binding_tool == "auto"
        assert inp.llm_tier == "none"
        assert inp.formats == ["tsv", "json"]

    def test_rank_input_accepts_custom_values(self):
        from dogneo.app.rank_pipeline import RankInput

        inp = RankInput(
            vcf_path=Path("my.vcf"),
            sample_id="DOG_001",
            alleles=["DLA-88*001:01"],
            binding_tool="iedb",
            formats=["tsv", "json", "fasta"],
        )
        assert inp.sample_id == "DOG_001"
        assert inp.alleles == ["DLA-88*001:01"]
        assert inp.binding_tool == "iedb"


class TestRankResult:
    """Test RankResult dataclass."""

    def test_rank_result_has_required_fields(self):
        from dogneo.app.rank_pipeline import RankResult

        result = RankResult(
            candidates=[],
            variants_total=8,
            variants_coding=6,
            peptides_total=228,
            binding_tool_used="iedb",
            explanations={},
        )
        assert result.variants_total == 8
        assert result.variants_coding == 6
        assert result.peptides_total == 228
        assert result.binding_tool_used == "iedb"
        assert result.candidates == []
        assert result.explanations == {}


class TestRunRankPipeline:
    """Test the main run_rank_pipeline function."""

    def test_pipeline_returns_rank_result(self, tmp_path):
        """Pipeline should return a RankResult with candidates."""
        # Use bundled demo VCF
        import dogneo.data as _data_pkg
        from dogneo.app.rank_pipeline import RankInput, RankResult, run_rank_pipeline
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST",
            binding_tool="none",
            formats=["tsv"],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = run_rank_pipeline(inp, output_dir)

        assert isinstance(result, RankResult)
        assert result.variants_total > 0
        assert result.variants_coding > 0
        assert result.binding_tool_used == "none"

    def test_pipeline_loads_bundled_alleles_when_none_provided(self, tmp_path):
        """When no alleles specified, should auto-load bundled DLA alleles."""
        import dogneo.data as _data_pkg
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST",
            alleles=[],  # empty — should auto-load
            binding_tool="none",
            formats=[],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = run_rank_pipeline(inp, output_dir)

        # Should have auto-loaded alleles (result tracks them)
        assert result.alleles_used is not None
        assert len(result.alleles_used) > 0

    def test_pipeline_exports_tsv(self, tmp_path):
        """Pipeline should create TSV output when requested."""
        import dogneo.data as _data_pkg
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_EXPORT",
            binding_tool="none",
            formats=["tsv"],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        run_rank_pipeline(inp, output_dir)

        tsv_path = output_dir / "TEST_EXPORT" / "candidates.tsv"
        assert tsv_path.exists()
        content = tsv_path.read_text()
        assert "gene" in content
        assert "mutant_peptide" in content

    def test_pipeline_exports_json(self, tmp_path):
        """Pipeline should create JSON output when requested."""
        import json

        import dogneo.data as _data_pkg
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_JSON",
            binding_tool="none",
            formats=["json"],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        run_rank_pipeline(inp, output_dir)

        json_path = output_dir / "TEST_JSON" / "candidates.json"
        assert json_path.exists()
        data = json.loads(json_path.read_text())
        assert "metadata" in data
        assert "candidates" in data
        assert data["metadata"]["sample_id"] == "TEST_JSON"

    def test_pipeline_uses_cached_proteome(self, tmp_path):
        """Pipeline should auto-detect cached proteome without re-download."""
        import dogneo.data as _data_pkg
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline
        from dogneo.data.manager import ReferenceDataManager
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        mgr = ReferenceDataManager()
        cached = mgr.get_proteome_path()
        if not cached:
            pytest.skip("Proteome not cached — run dogneo setup first")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_CACHE",
            binding_tool="none",
            formats=[],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        result = run_rank_pipeline(inp, output_dir)

        # Should have generated peptides (proves proteome was loaded)
        assert result.peptides_total > 0
