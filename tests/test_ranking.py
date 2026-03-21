"""Tests for dogneo.core.ranking — scoring functions and candidate ranking."""
from __future__ import annotations

import pytest

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.ranking import (
    NeoantigenCandidate,
    ScoringWeights,
    _score_binding,
    _score_expression,
    _score_vaf,
    _score_self_difference,
    _score_agretopicity,
    _score_caller_agreement,
    rank_candidates,
    build_candidates,
)
from dogneo.core.variants import SomaticVariant


# ---------------------------------------------------------------------------
# Individual scoring functions
# ---------------------------------------------------------------------------

class TestScoreBinding:
    """Tests for binding affinity scoring."""

    def test_strong_binder(self):
        """< 50 nM → score ≈ 1.0."""
        score = _score_binding(30.0)
        assert score > 0.9

    def test_weak_binder(self):
        """> 5000 nM → score ≈ 0.0."""
        score = _score_binding(6000.0)
        assert score < 0.1

    def test_medium_binder(self):
        """500 nM → score ≈ 0.5."""
        score = _score_binding(500.0)
        assert 0.3 < score < 0.7

    def test_zero_affinity(self):
        """Edge case: 0 nM → score = 1.0."""
        assert _score_binding(0.0) == 1.0

    def test_negative_affinity(self):
        """Edge case: negative → score = 1.0."""
        assert _score_binding(-10.0) == 1.0


class TestScoreExpression:
    """Tests for expression level scoring."""

    def test_high_tpm(self):
        """TPM 100 → score ≈ 1.0."""
        score = _score_expression(100.0)
        assert score > 0.9

    def test_low_tpm(self):
        """TPM 0.5 → score near 0."""
        score = _score_expression(0.5)
        assert score < 0.2

    def test_zero_tpm(self):
        assert _score_expression(0.0) == 0.0

    def test_negative_tpm(self):
        assert _score_expression(-1.0) == 0.0


class TestScoreVaf:
    """Tests for variant allele frequency scoring."""

    def test_high_vaf(self):
        """VAF 0.5 → capped at 1.0."""
        assert _score_vaf(0.5) == 1.0

    def test_above_half(self):
        """VAF > 0.5 → still 1.0."""
        assert _score_vaf(0.8) == 1.0

    def test_low_vaf(self):
        """VAF 0.1 → 0.2."""
        assert _score_vaf(0.1) == pytest.approx(0.2)

    def test_zero_vaf(self):
        assert _score_vaf(0.0) == 0.0


class TestScoreSelfDifference:
    """Tests for wild-type vs mutant difference scoring."""

    def test_single_aa_change(self):
        """One amino acid difference out of 9 → 1/3 ≈ 0.333."""
        score = _score_self_difference("AAAAAAAAA", "AAAAAAAAB")
        assert score == pytest.approx(1 / 3, abs=0.01)

    def test_length_mismatch(self):
        """Different lengths → max score 1.0."""
        score = _score_self_difference("AAAA", "AAAAA")
        assert score == 1.0

    def test_identical(self):
        score = _score_self_difference("AAAA", "AAAA")
        assert score == 0.0


class TestScoreAgretopicity:
    """Tests for agretopicity (WT/mut binding ratio) scoring."""

    def test_uses_wt_binding(self):
        """wt=500, mut=50 → ratio=10 → high agretopicity."""
        score = _score_agretopicity(wt_binding_nm=500.0, mut_binding_nm=50.0)
        assert score > 0.8

    def test_neutral_when_equal(self):
        """wt=100, mut=100 → ratio=1 → score ≈ 0.5."""
        score = _score_agretopicity(wt_binding_nm=100.0, mut_binding_nm=100.0)
        assert score == pytest.approx(0.5, abs=0.05)

    def test_low_agretopicity(self):
        """wt=10, mut=1000 → ratio=0.01 → low score."""
        score = _score_agretopicity(wt_binding_nm=10.0, mut_binding_nm=1000.0)
        assert score < 0.3

    def test_zero_mut_binding(self):
        """mut=0 → score = 1.0."""
        assert _score_agretopicity(wt_binding_nm=100.0, mut_binding_nm=0.0) == 1.0

    def test_zero_wt_binding(self):
        """wt=0 → assume very weak WT → high agretopicity."""
        score = _score_agretopicity(wt_binding_nm=0.0, mut_binding_nm=50.0)
        assert score > 0.8


class TestScoreCallerAgreement:
    """Tests for caller agreement scoring."""

    def test_single_caller(self):
        assert _score_caller_agreement(1) == pytest.approx(1 / 3, abs=0.01)

    def test_all_callers(self):
        assert _score_caller_agreement(3) == 1.0

    def test_zero_callers(self):
        assert _score_caller_agreement(0) == 0.0


# ---------------------------------------------------------------------------
# NeoantigenCandidate
# ---------------------------------------------------------------------------

class TestNeoantigenCandidate:
    """Tests for NeoantigenCandidate dataclass."""

    def test_candidate_id(self, sample_candidate: NeoantigenCandidate):
        cid = sample_candidate.candidate_id
        assert "chr1:15000100:G>A" in cid
        assert "RAIVGAPPS" in cid
        assert "DLA-88*001:01" in cid

    def test_to_dict(self, sample_candidate: NeoantigenCandidate):
        d = sample_candidate.to_dict()
        assert d["gene"] == "TP53"
        assert d["rank"] == 1
        assert d["mutant_peptide"] == "RAIVGAPPS"
        assert d["allele"] == "DLA-88*001:01"
        assert d["binding_affinity_nm"] == pytest.approx(42.5)
        assert "composite_score" in d
        assert "score_components" in d


# ---------------------------------------------------------------------------
# rank_candidates
# ---------------------------------------------------------------------------

def _make_candidate(affinity_nm: float, vaf: float, tpm: float) -> NeoantigenCandidate:
    """Helper to build a candidate with specific values for ranking tests."""
    variant = SomaticVariant(
        chrom="chr1", pos=100, ref="A", alt="T",
        gene="TEST", effect="missense_variant",
        vaf=vaf, depth=50, expression_tpm=tpm,
    )
    peptide = MutantPeptide(
        gene="TEST", variant_id=variant.variant_id,
        mutation="p.A100T", wt_sequence="AAAAAAAAA",
        mut_sequence="AAAAAAAAB", position=8, length=9,
    )
    binding = BindingPrediction(
        peptide_sequence="AAAAAAAAB", allele="DLA-88*001:01",
        affinity_nm=affinity_nm, percentile_rank=5.0,
    )
    return NeoantigenCandidate(
        variant=variant, peptide=peptide, binding=binding,
        expression_tpm=tpm,
    )


class TestRankCandidates:
    """Tests for the rank_candidates function."""

    def test_rank_order(self):
        """Candidate with stronger binding + higher expression should rank first."""
        strong = _make_candidate(affinity_nm=30.0, vaf=0.5, tpm=100.0)
        weak = _make_candidate(affinity_nm=3000.0, vaf=0.1, tpm=1.0)
        ranked = rank_candidates([weak, strong])
        assert ranked[0].rank == 1
        assert ranked[0].composite_score > ranked[1].composite_score

    def test_scores_populated(self):
        c = _make_candidate(affinity_nm=100.0, vaf=0.3, tpm=50.0)
        ranked = rank_candidates([c])
        assert ranked[0].composite_score > 0
        assert "binding" in ranked[0].score_components

    def test_empty_list(self):
        result = rank_candidates([])
        assert result == []


# ---------------------------------------------------------------------------
# build_candidates
# ---------------------------------------------------------------------------

class TestBuildCandidates:
    """Tests for assembling candidates from components."""

    def test_build_matches(
        self,
        sample_variant: SomaticVariant,
        sample_peptide: MutantPeptide,
        sample_binding: BindingPrediction,
    ):
        peptides_map = {sample_variant.variant_id: [sample_peptide]}
        preds_map = {sample_peptide.mut_sequence: [sample_binding]}
        result = build_candidates([sample_variant], peptides_map, preds_map)
        assert len(result) == 1
        assert result[0].variant.gene == "TP53"
        assert result[0].binding.allele == "DLA-88*001:01"

    def test_build_empty(self):
        result = build_candidates([], {}, {})
        assert result == []

    def test_build_no_matching_predictions(
        self,
        sample_variant: SomaticVariant,
        sample_peptide: MutantPeptide,
    ):
        peptides_map = {sample_variant.variant_id: [sample_peptide]}
        preds_map: dict = {}  # no predictions
        result = build_candidates([sample_variant], peptides_map, preds_map)
        assert len(result) == 0
