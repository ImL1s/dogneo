"""Tests for the LLM pipeline explainer.

TDD Phase 4: Step-by-step plain-language explanations at each pipeline stage.
"""
from __future__ import annotations

from unittest.mock import MagicMock

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.core.variants import SomaticVariant


def _make_variant(gene: str = "TP53", vaf: float = 0.4, tpm: float = 85.0) -> SomaticVariant:
    return SomaticVariant(
        chrom="chr1", pos=100, ref="G", alt="A",
        gene=gene, effect="missense_variant", hgvs_p="p.V245I",
        vaf=vaf, depth=50, alt_depth=20, filter_status="PASS",
        expression_tpm=tpm,
    )


def _make_candidate(
    gene: str = "TP53", affinity: float = 42.5, tpm: float = 85.0
) -> NeoantigenCandidate:
    variant = _make_variant(gene=gene, tpm=tpm)
    peptide = MutantPeptide(
        gene=gene, variant_id="chr1:100:G>A", mutation="p.V245I",
        wt_sequence="RAVVGAPPS", mut_sequence="RAIVGAPPS",
        position=2, length=9, mhc_class=1,
    )
    binding = BindingPrediction(
        peptide_sequence="RAIVGAPPS", allele="DLA-88*001:01",
        affinity_nm=affinity, percentile_rank=0.8,
        tool="iedb", mhc_class=1,
    )
    return NeoantigenCandidate(
        variant=variant, peptide=peptide, binding=binding,
        expression_tpm=tpm, composite_score=0.85, rank=1,
    )


class TestPipelineExplainer:

    def test_explainer_can_be_instantiated(self):
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        explainer = PipelineExplainer(router=mock_router)
        assert explainer is not None

    def test_explain_variants_mentions_gene_names(self):
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        mock_router.generate.return_value = "Found mutations in TP53 and BRAF genes."

        explainer = PipelineExplainer(router=mock_router)
        variants = [_make_variant("TP53"), _make_variant("BRAF")]
        coding = variants  # all are coding

        result = explainer.explain_variants(variants, coding)
        assert isinstance(result, str)
        assert len(result) > 0
        # LLM was called
        mock_router.generate.assert_called_once()

    def test_explain_binding_mentions_strong_binders(self):
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        mock_router.generate.return_value = "3 peptides are strong binders (<50nM)."

        explainer = PipelineExplainer(router=mock_router)
        candidates = [
            _make_candidate("TP53", affinity=42.5),
            _make_candidate("BRAF", affinity=350.0),
        ]

        result = explainer.explain_binding(candidates)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_explain_ranking_describes_top_candidate(self):
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        mock_router.generate.return_value = "Top candidate is TP53 V245I with strong binding."

        explainer = PipelineExplainer(router=mock_router)
        candidates = [_make_candidate("TP53"), _make_candidate("BRAF", affinity=350.0)]

        result = explainer.explain_ranking(candidates[:5])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_explain_for_owner_is_plain_language(self):
        from dogneo.app.rank_pipeline import RankResult
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        mock_router.generate.return_value = (
            "We analyzed your dog's tumor DNA and found 6 mutations. "
            "The top candidate targets the TP53 gene."
        )

        explainer = PipelineExplainer(router=mock_router)
        rank_result = RankResult(
            candidates=[_make_candidate()],
            variants_total=8, variants_coding=6,
            peptides_total=228, binding_tool_used="iedb",
        )

        result = explainer.explain_for_owner(rank_result)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_explainer_graceful_when_llm_unavailable(self):
        from dogneo.llm.explainer import PipelineExplainer

        mock_router = MagicMock()
        mock_router.generate.side_effect = RuntimeError("All backends failed")

        explainer = PipelineExplainer(router=mock_router)
        variants = [_make_variant("TP53")]

        # Should not raise, should return fallback message
        result = explainer.explain_variants(variants, variants)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_explainer_works_without_router(self):
        """When no LLM router provided, return structured summary without AI."""
        from dogneo.llm.explainer import PipelineExplainer

        explainer = PipelineExplainer(router=None)
        variants = [_make_variant("TP53"), _make_variant("BRAF")]

        result = explainer.explain_variants(variants, variants)
        assert isinstance(result, str)
        assert "TP53" in result
        assert "BRAF" in result
