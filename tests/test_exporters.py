"""Tests for dogneo.export.exporters — FASTA, TSV, JSON output."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from dogneo import __version__, RUO_DISCLAIMER
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.export.exporters import export_fasta, export_tsv, export_json


# ---------------------------------------------------------------------------
# FASTA export
# ---------------------------------------------------------------------------

class TestExportFasta:
    """Tests for FASTA peptide export."""

    def test_fasta_header_content(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        """FASTA header should contain gene, mutation, allele, rank, score."""
        out = tmp_output_dir / "test.fasta"
        export_fasta([sample_candidate], out)

        text = out.read_text()
        assert ">TP53_p.Val245Ile" in text
        assert "allele=DLA-88*001:01" in text
        assert "rank=1" in text

    def test_fasta_ruo_disclaimer(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        """FASTA should contain RUO disclaimer."""
        out = tmp_output_dir / "test.fasta"
        export_fasta([sample_candidate], out)
        text = out.read_text()
        assert "RESEARCH USE ONLY" in text

    def test_fasta_sequence(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        """FASTA should contain the mutant peptide sequence."""
        out = tmp_output_dir / "test.fasta"
        export_fasta([sample_candidate], out)
        text = out.read_text()
        assert "RAIVGAPPS" in text

    def test_fasta_top_n(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        """top_n=0 → should write all (no limit)."""
        out = tmp_output_dir / "top.fasta"
        export_fasta([sample_candidate], out, top_n=1)
        lines = [l for l in out.read_text().splitlines() if l.startswith(">")]
        assert len(lines) == 1

    def test_fasta_returns_path(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "ret.fasta"
        result = export_fasta([sample_candidate], out)
        assert result == out


# ---------------------------------------------------------------------------
# TSV export
# ---------------------------------------------------------------------------

class TestExportTsv:
    """Tests for TSV tabular export."""

    def test_tsv_columns(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        """TSV should have expected column headers."""
        out = tmp_output_dir / "test.tsv"
        export_tsv([sample_candidate], out)

        # Read non-comment lines
        with open(out) as f:
            lines = [l for l in f if not l.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        row = next(reader)

        expected_cols = {
            "rank", "gene", "mutation", "mutant_peptide",
            "allele", "binding_affinity_nm", "composite_score",
        }
        assert expected_cols.issubset(set(row.keys()))

    def test_tsv_ruo_comment(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "test.tsv"
        export_tsv([sample_candidate], out)
        first_line = out.read_text().splitlines()[0]
        assert first_line.startswith("# DogNeo")
        assert "RESEARCH USE ONLY" in first_line

    def test_tsv_data_values(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "vals.tsv"
        export_tsv([sample_candidate], out)
        with open(out) as f:
            lines = [l for l in f if not l.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        row = next(reader)
        assert row["gene"] == "TP53"
        assert row["mutant_peptide"] == "RAIVGAPPS"


# ---------------------------------------------------------------------------
# JSON export
# ---------------------------------------------------------------------------

class TestExportJson:
    """Tests for JSON structured export."""

    def test_json_metadata(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "test.json"
        export_json([sample_candidate], out, sample_id="DOG001")

        data = json.loads(out.read_text())
        assert "metadata" in data
        assert "candidates" in data
        assert data["metadata"]["tool"] == "DogNeo"
        assert data["metadata"]["version"] == __version__
        assert data["metadata"]["sample_id"] == "DOG001"
        assert "RESEARCH USE ONLY" in data["metadata"]["disclaimer"]

    def test_json_candidate_content(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "test.json"
        export_json([sample_candidate], out)

        data = json.loads(out.read_text())
        c = data["candidates"][0]
        assert c["gene"] == "TP53"
        assert c["mutant_peptide"] == "RAIVGAPPS"
        assert c["allele"] == "DLA-88*001:01"

    def test_json_top_n(
        self,
        sample_candidate: NeoantigenCandidate,
        tmp_output_dir: Path,
    ):
        out = tmp_output_dir / "topn.json"
        export_json([sample_candidate, sample_candidate], out, top_n=1)
        data = json.loads(out.read_text())
        assert data["metadata"]["exported_candidates"] == 1
        assert data["metadata"]["total_candidates"] == 2
        assert len(data["candidates"]) == 1
