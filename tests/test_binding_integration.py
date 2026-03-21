"""Tests for IEDB integration into rank pipeline.

TDD Phase 2.2: rank pipeline uses IEDB when binding_tool is auto/iedb.
"""
from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


MOCK_IEDB_RESPONSE = (
    "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tpercentile_rank\tic50\n"
    "DLA-88*001:01\t1\t1\t9\t9\tRAIVGAPPS\tnetmhcpan_el\t1.5\t125.0\n"
    "DLA-88*001:01\t2\t1\t9\t9\tYLVSGAPPS\tnetmhcpan_el\t0.3\t42.0\n"
)


class TestPipelineWithIEDB:
    """Test rank pipeline using IEDB binding predictions."""

    def _make_mock_response(self, peptides: list[str], alleles: list[str]) -> MagicMock:
        """Create a mock IEDB API response with predictions for given peptides."""
        lines = ["allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tpercentile_rank\tic50\n"]
        for i, pep in enumerate(peptides):
            for allele in alleles:
                affinity = 50.0 + (i * 100)  # varied affinities
                rank = 0.5 + (i * 0.5)
                lines.append(
                    f"{allele}\t{i+1}\t1\t{len(pep)}\t{len(pep)}\t{pep}\tnetmhcpan_el\t{rank}\t{affinity}\n"
                )
        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = "".join(lines)
        return mock_resp

    def test_pipeline_with_iedb_produces_scored_candidates(self, tmp_path):
        """When binding_tool=iedb, candidates should have real affinity values."""
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline

        import dogneo.data as _data_pkg
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_IEDB",
            binding_tool="iedb",
            formats=[],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Dynamic mock: intercept the actual POST, extract peptides, return matching predictions
        def dynamic_iedb_response(*args, **kwargs):
            data = kwargs.get("data", {})
            peptides = data.get("sequence_text", "").strip().split("\n")
            alleles_list = data.get("allele", "").split(",")
            return self._make_mock_response(peptides, alleles_list)

        # Use isolated cache dir to avoid stale cache from previous runs
        with patch("dogneo.core.iedb_client.DEFAULT_CACHE_DIR", tmp_path / "iedb_cache"), \
             patch("dogneo.core.iedb_client.requests.post", side_effect=dynamic_iedb_response):
            result = run_rank_pipeline(inp, output_dir)

        assert result.binding_tool_used == "iedb"
        # At least some candidates should have real affinity (not NaN)
        scored = [c for c in result.candidates if not math.isnan(c.binding.affinity_nm)]
        assert len(scored) > 0

    def test_pipeline_with_iedb_ranks_by_score(self, tmp_path):
        """Scored candidates should be ranked by composite_score descending."""
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline

        import dogneo.data as _data_pkg
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_RANKED",
            binding_tool="iedb",
            formats=[],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        def dynamic_iedb_response(*args, **kwargs):
            data = kwargs.get("data", {})
            peptides = data.get("sequence_text", "").strip().split("\n")
            alleles_list = data.get("allele", "").split(",")
            return self._make_mock_response(peptides, alleles_list)

        with patch("dogneo.core.iedb_client.DEFAULT_CACHE_DIR", tmp_path / "iedb_cache"), \
             patch("dogneo.core.iedb_client.requests.post", side_effect=dynamic_iedb_response):
            result = run_rank_pipeline(inp, output_dir)

        scored = [c for c in result.candidates if not math.isnan(c.binding.affinity_nm)]
        assert len(scored) >= 2
        # Should be sorted descending by composite_score
        scores = [c.composite_score for c in scored]
        assert scores == sorted(scores, reverse=True)

    def test_pipeline_auto_resolves_to_iedb(self, tmp_path):
        """binding_tool=auto should resolve to iedb when netMHCpan not installed."""
        from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline

        import dogneo.data as _data_pkg
        demo_dir = Path(_data_pkg.__file__).parent / "demo"
        vcf_path = demo_dir / "canine_osteosarcoma.vcf"

        if not vcf_path.exists():
            pytest.skip("Demo VCF not available")

        inp = RankInput(
            vcf_path=vcf_path,
            sample_id="TEST_AUTO",
            binding_tool="auto",
            formats=[],
        )
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("shutil.which", return_value=None), \
             patch("dogneo.core.iedb_client.requests.post") as mock_post:
            mock_post.return_value = self._make_mock_response(
                ["RAIVGAPPS"], ["DLA-88*001:01"]
            )
            result = run_rank_pipeline(inp, output_dir)

        assert result.binding_tool_used == "iedb"
