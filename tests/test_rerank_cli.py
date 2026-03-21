"""Tests for dogneo rerank CLI command.

TDD Phase 3: CLI integration for rerank.
"""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from dogneo.cli import cli


def _create_test_files(tmp_path: Path) -> tuple[Path, Path]:
    """Create candidates.json and binding.tsv for CLI testing."""
    candidates = {
        "metadata": {
            "tool": "DogNeo", "version": "0.1.0",
            "sample_id": "TEST", "total_candidates": 1, "exported_candidates": 1,
        },
        "candidates": [{
            "rank": 1, "gene": "TP53", "variant_id": "chr1:100:G>A",
            "mutation": "p.V245I", "mutant_peptide": "RAIVGAPPS",
            "wildtype_peptide": "RAVVGAPPS", "peptide_length": 9,
            "mhc_class": 1, "allele": "DLA-88*001:01",
            "binding_affinity_nm": float("nan"), "binding_percentile": float("nan"),
            "binding_tool": "pending", "vaf": 0.4, "expression_tpm": 85.3,
            "composite_score": 0.0, "score_components": {}, "num_callers": 0,
        }],
    }
    candidates_path = tmp_path / "candidates.json"
    candidates_path.write_text(json.dumps(candidates))

    binding_path = tmp_path / "binding.tsv"
    binding_path.write_text(
        "peptide\tallele\taffinity_nm\tpercentile_rank\n"
        "RAIVGAPPS\tDLA-88*001:01\t42.5\t0.8\n"
    )
    return candidates_path, binding_path


class TestRerankCLI:

    def test_rerank_command_exists(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["rerank", "--help"])
        assert result.exit_code == 0
        assert "candidates" in result.output.lower()

    def test_rerank_produces_output(self, tmp_path):
        candidates_path, binding_path = _create_test_files(tmp_path)
        output_dir = tmp_path / "reranked"

        runner = CliRunner()
        result = runner.invoke(cli, [
            "rerank",
            "--candidates", str(candidates_path),
            "--binding", str(binding_path),
            "--output-dir", str(output_dir),
        ])

        assert result.exit_code == 0
        assert "Re-ranked" in result.output or "reranked" in result.output.lower()
        assert (output_dir / "candidates.tsv").exists()
