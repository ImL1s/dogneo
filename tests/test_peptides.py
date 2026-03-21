"""Tests for dogneo.core.peptides — peptide generation and protein DB."""
from __future__ import annotations

import pytest

from dogneo.core.peptides import MutantPeptide, ProteinDatabase


# ---------------------------------------------------------------------------
# MutantPeptide
# ---------------------------------------------------------------------------

class TestMutantPeptide:
    """Tests for MutantPeptide dataclass and methods."""

    @pytest.mark.parametrize("length", [8, 9, 10, 11])
    def test_peptide_length(self, length: int):
        """Peptide object should store the declared length."""
        pep = MutantPeptide(
            gene="TP53",
            variant_id="chr1:100:A>T",
            mutation="p.V245I",
            wt_sequence="A" * length,
            mut_sequence="A" * (length - 1) + "T",
            position=length - 1,
            length=length,
        )
        assert pep.length == length
        assert len(pep.mut_sequence) == length

    def test_peptide_id_format(self):
        pep = MutantPeptide(
            gene="BRAF",
            variant_id="chr5:200:T>C",
            mutation="p.V600A",
            wt_sequence="AAAAAAAAV",
            mut_sequence="AAAAAAAAA",
            position=8,
            length=9,
        )
        assert pep.peptide_id == "chr5:200:T>C|AAAAAAAAA|9mer"

    def test_self_different(self, sample_peptide: MutantPeptide):
        """wt=RAVVGAPPS, mut=RAIVGAPPS → different at position 2."""
        assert sample_peptide.is_self_different() is True

    def test_self_identical(self):
        pep = MutantPeptide(
            gene="X",
            variant_id="chr1:1:A>A",
            mutation="p.=",
            wt_sequence="AAAA",
            mut_sequence="AAAA",
            position=0,
            length=4,
        )
        assert pep.is_self_different() is False

    def test_mutation_position_in_window(self, sample_peptide: MutantPeptide):
        """Mutation should be at the declared position within the peptide window."""
        pos = sample_peptide.position
        assert 0 <= pos < sample_peptide.length
        assert sample_peptide.wt_sequence[pos] != sample_peptide.mut_sequence[pos]


# ---------------------------------------------------------------------------
# ProteinDatabase
# ---------------------------------------------------------------------------

class TestProteinDatabase:
    """Tests for ProteinDatabase FASTA loading."""

    def test_load_fasta(self, tmp_path):
        """ProteinDatabase should load a minimal protein FASTA."""
        fasta = tmp_path / "proteins.fa"
        fasta.write_text(
            ">ENSCAFT00000013077 gene=TP53\n"
            "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
            ">ENSCAFT00000005678 gene=BRAF\n"
            "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWN\n"
        )
        db = ProteinDatabase()
        db.load_fasta(fasta)
        assert db.get_sequence(gene="TP53") is not None
        assert db.get_sequence(transcript_id="ENSCAFT00000005678") is not None

    def test_empty_db_returns_none(self):
        db = ProteinDatabase()
        assert db.get_sequence(gene="NONEXISTENT") is None

    def test_get_by_transcript_id(self, tmp_path):
        fasta = tmp_path / "p.fa"
        fasta.write_text(">TX001 gene=GENE1\nMAAA\n")
        db = ProteinDatabase()
        db.load_fasta(fasta)
        seq = db.get_sequence(transcript_id="TX001")
        assert seq == "MAAA"
