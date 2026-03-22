"""Tests for the rerank pipeline.

TDD Phase 3: Import external binding results and re-score candidates.
"""
from __future__ import annotations

import json
from pathlib import Path


def _create_candidates_json(path: Path, sample_id: str = "TEST") -> Path:
    """Create a minimal candidates.json for testing."""
    data = {
        "metadata": {
            "tool": "DogNeo",
            "version": "0.1.0",
            "sample_id": sample_id,
            "total_candidates": 2,
            "exported_candidates": 2,
        },
        "candidates": [
            {
                "rank": 1,
                "gene": "TP53",
                "variant_id": "chr1:15000100:G>A",
                "mutation": "p.Val245Ile",
                "mutant_peptide": "RAIVGAPPS",
                "wildtype_peptide": "RAVVGAPPS",
                "peptide_length": 9,
                "mhc_class": 1,
                "allele": "DLA-88*001:01",
                "binding_affinity_nm": float("nan"),
                "binding_percentile": float("nan"),
                "binding_tool": "pending",
                "vaf": 0.40,
                "expression_tpm": 85.3,
                "composite_score": 0.0,
                "score_components": {},
                "num_callers": 0,
            },
            {
                "rank": 2,
                "gene": "BRAF",
                "variant_id": "chr5:32500200:T>C",
                "mutation": "p.Val600Ala",
                "mutant_peptide": "YLVSGAPPS",
                "wildtype_peptide": "YLVSGATPS",
                "peptide_length": 9,
                "mhc_class": 1,
                "allele": "DLA-88*001:01",
                "binding_affinity_nm": float("nan"),
                "binding_percentile": float("nan"),
                "binding_tool": "pending",
                "vaf": 0.375,
                "expression_tpm": 2.1,
                "composite_score": 0.0,
                "score_components": {},
                "num_callers": 0,
            },
        ],
    }
    path.write_text(json.dumps(data))
    return path


def _create_binding_tsv(path: Path) -> Path:
    """Create a generic binding results TSV."""
    lines = [
        "peptide\tallele\taffinity_nm\tpercentile_rank",
        "RAIVGAPPS\tDLA-88*001:01\t42.5\t0.8",
        "YLVSGAPPS\tDLA-88*001:01\t350.0\t3.2",
    ]
    path.write_text("\n".join(lines))
    return path


def _create_netmhcpan_tsv(path: Path) -> Path:
    """Create a NetMHCpan-style XLS/TSV output."""
    lines = [
        "Pos\tMHC\tPeptide\tCore\tOf\tGp\tGl\tIp\tIl\tIcore\tIdentity\tScore_EL\t%Rank_EL\tScore_BA\tAff(nM)\t%Rank_BA\tBindLevel",
        "1\tDLA-88*001:01\tRAIVGAPPS\tRAIVGAPPS\t0\t0\t0\t0\t0\tRAIVGAPPS\tTP53\t0.85\t0.5\t0.92\t42.5\t0.8\t<=SB",
        "1\tDLA-88*001:01\tYLVSGAPPS\tYLVSGAPPS\t0\t0\t0\t0\t0\tYLVSGAPPS\tBRAF\t0.45\t2.1\t0.55\t350.0\t3.2\t<=WB",
    ]
    path.write_text("\n".join(lines))
    return path


class TestRerankInput:
    """Test RerankInput dataclass."""

    def test_rerank_input_has_required_fields(self):
        from dogneo.app.rerank_pipeline import RerankInput

        inp = RerankInput(
            candidates_path=Path("candidates.json"),
            binding_path=Path("binding.tsv"),
        )
        assert inp.candidates_path == Path("candidates.json")
        assert inp.binding_path == Path("binding.tsv")
        assert inp.binding_format == "auto"


class TestBindingFormatDetection:
    """Test automatic detection of binding result formats."""

    def test_detects_generic_tsv(self, tmp_path):
        from dogneo.app.rerank_pipeline import detect_binding_format

        tsv_path = _create_binding_tsv(tmp_path / "binding.tsv")
        assert detect_binding_format(tsv_path) == "tsv"

    def test_detects_netmhcpan(self, tmp_path):
        from dogneo.app.rerank_pipeline import detect_binding_format

        tsv_path = _create_netmhcpan_tsv(tmp_path / "netmhcpan.tsv")
        assert detect_binding_format(tsv_path) == "netmhcpan"


class TestParseBindingResults:
    """Test parsing of different binding result formats."""

    def test_parse_generic_tsv(self, tmp_path):
        from dogneo.app.rerank_pipeline import parse_binding_results

        tsv_path = _create_binding_tsv(tmp_path / "binding.tsv")
        results = parse_binding_results(tsv_path, fmt="tsv")

        assert len(results) == 2
        assert ("RAIVGAPPS", "DLA-88*001:01") in results
        assert results[("RAIVGAPPS", "DLA-88*001:01")].affinity_nm == 42.5

    def test_parse_netmhcpan(self, tmp_path):
        from dogneo.app.rerank_pipeline import parse_binding_results

        tsv_path = _create_netmhcpan_tsv(tmp_path / "netmhcpan.tsv")
        results = parse_binding_results(tsv_path, fmt="netmhcpan")

        assert len(results) == 2
        assert ("RAIVGAPPS", "DLA-88*001:01") in results
        pred = results[("RAIVGAPPS", "DLA-88*001:01")]
        assert pred.affinity_nm == 42.5
        assert pred.tool == "netmhcpan"


class TestRunRerankPipeline:
    """Test the full rerank pipeline."""

    def test_rerank_merges_binding_and_rescores(self, tmp_path):
        from dogneo.app.rerank_pipeline import RerankInput, run_rerank_pipeline

        candidates_path = _create_candidates_json(tmp_path / "candidates.json")
        binding_path = _create_binding_tsv(tmp_path / "binding.tsv")
        output_dir = tmp_path / "reranked"
        output_dir.mkdir()

        inp = RerankInput(
            candidates_path=candidates_path,
            binding_path=binding_path,
        )
        result = run_rerank_pipeline(inp, output_dir)

        # Should have re-scored candidates
        assert len(result.candidates) == 2
        # TP53 with 42.5 nM should rank higher than BRAF with 350 nM
        assert result.candidates[0].binding.affinity_nm == 42.5
        assert result.candidates[0].composite_score > result.candidates[1].composite_score

    def test_rerank_exports_results(self, tmp_path):
        from dogneo.app.rerank_pipeline import RerankInput, run_rerank_pipeline

        candidates_path = _create_candidates_json(tmp_path / "candidates.json")
        binding_path = _create_binding_tsv(tmp_path / "binding.tsv")
        output_dir = tmp_path / "reranked"
        output_dir.mkdir()

        inp = RerankInput(
            candidates_path=candidates_path,
            binding_path=binding_path,
            formats=["tsv", "json"],
        )
        result = run_rerank_pipeline(inp, output_dir)

        assert (output_dir / "candidates.tsv").exists()
        assert (output_dir / "candidates.json").exists()
