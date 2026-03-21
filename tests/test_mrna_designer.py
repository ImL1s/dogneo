"""Tests for the mRNA construct designer.

TDD Phase 5: Codon-optimized mRNA construct from top candidates.
"""
from __future__ import annotations

import pytest

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.core.variants import SomaticVariant


def _make_candidate(gene: str, peptide_seq: str, rank: int = 1) -> NeoantigenCandidate:
    return NeoantigenCandidate(
        variant=SomaticVariant(
            chrom="chr1", pos=100, ref="G", alt="A",
            gene=gene, effect="missense_variant", hgvs_p="p.V245I",
            vaf=0.4, depth=50, alt_depth=20, filter_status="PASS",
            expression_tpm=85.0,
        ),
        peptide=MutantPeptide(
            gene=gene, variant_id=f"chr1:100:G>A",
            mutation="p.V245I", wt_sequence="A" * len(peptide_seq),
            mut_sequence=peptide_seq, position=2,
            length=len(peptide_seq), mhc_class=1,
        ),
        binding=BindingPrediction(
            peptide_sequence=peptide_seq, allele="DLA-88*001:01",
            affinity_nm=42.5, percentile_rank=0.8,
            tool="iedb", mhc_class=1,
        ),
        expression_tpm=85.0, composite_score=0.85, rank=rank,
    )


class TestCodonOptimize:

    def test_codon_optimize_returns_dna_sequence(self):
        from dogneo.core.mrna_designer import codon_optimize

        result = codon_optimize("MKVL")
        assert isinstance(result, str)
        assert len(result) == 4 * 3  # 4 amino acids × 3 nucleotides
        assert all(c in "ACGU" for c in result)

    def test_codon_optimize_output_translates_back(self):
        from dogneo.core.mrna_designer import codon_optimize, _translate_rna

        protein = "RAIVGAPPS"
        rna = codon_optimize(protein)
        translated = _translate_rna(rna)
        assert translated == protein

    def test_codon_optimize_uses_canine_table(self):
        from dogneo.core.mrna_designer import codon_optimize, CANINE_CODON_TABLE

        # Canine table should exist and have entries
        assert len(CANINE_CODON_TABLE) == 20  # 20 amino acids
        # Each entry should be a valid RNA codon
        for aa, codon in CANINE_CODON_TABLE.items():
            assert len(codon) == 3
            assert all(c in "ACGU" for c in codon)


class TestMRNAConstruct:

    def test_construct_has_utr_and_polya(self):
        from dogneo.core.mrna_designer import design_construct

        candidates = [_make_candidate("TP53", "RAIVGAPPS")]
        construct = design_construct(candidates, top_n=1)

        assert construct.five_prime_utr != ""
        assert construct.three_prime_utr != ""
        assert construct.poly_a_length >= 100

    def test_construct_includes_linkers_between_epitopes(self):
        from dogneo.core.mrna_designer import design_construct

        candidates = [
            _make_candidate("TP53", "RAIVGAPPS", rank=1),
            _make_candidate("BRAF", "YLVSGAPPS", rank=2),
        ]
        construct = design_construct(candidates, top_n=2)

        assert len(construct.linker_sequences) > 0
        assert construct.antigen_cassette != ""

    def test_construct_full_sequence_is_valid_rna(self):
        from dogneo.core.mrna_designer import design_construct

        candidates = [_make_candidate("TP53", "RAIVGAPPS")]
        construct = design_construct(candidates, top_n=1)

        seq = construct.full_sequence
        assert len(seq) > 0
        assert all(c in "ACGU" for c in seq)

    def test_construct_selects_top_n(self):
        from dogneo.core.mrna_designer import design_construct

        candidates = [
            _make_candidate("TP53", "RAIVGAPPS", rank=1),
            _make_candidate("BRAF", "YLVSGAPPS", rank=2),
            _make_candidate("KRAS", "MTEYKLVVV", rank=3),
        ]
        construct = design_construct(candidates, top_n=2)

        # Should only include top 2 epitopes
        assert construct.epitope_count == 2

    def test_construct_to_fasta(self):
        from dogneo.core.mrna_designer import design_construct

        candidates = [_make_candidate("TP53", "RAIVGAPPS")]
        construct = design_construct(candidates, top_n=1)

        fasta = construct.to_fasta()
        assert fasta.startswith(">")
        assert "\n" in fasta
        # Sequence part should be valid RNA
        seq_lines = [l for l in fasta.split("\n") if not l.startswith(">")]
        seq = "".join(seq_lines)
        assert all(c in "ACGU" for c in seq)
