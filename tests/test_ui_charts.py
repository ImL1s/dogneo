"""Tests for Plotly chart builders.

TDD Phase 6: Streamlit UI visualization components.
"""
from __future__ import annotations

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.core.variants import SomaticVariant


def _make_candidates(n: int = 5) -> list[NeoantigenCandidate]:
    candidates = []
    genes = ["TP53", "BRAF", "KRAS", "PIK3CA", "PTEN"]
    peptides = ["RAIVGAPPS", "YLVSGAPPS", "MTEYKLVVV", "HGDILHRAL", "FLTPKKLQC"]
    affinities = [42.5, 150.0, 350.0, 85.0, 1200.0]

    for i in range(min(n, 5)):
        candidates.append(NeoantigenCandidate(
            variant=SomaticVariant(
                chrom="chr1", pos=100 + i, ref="G", alt="A",
                gene=genes[i], effect="missense_variant", hgvs_p=f"p.V{i}I",
                vaf=0.4 - (i * 0.05), depth=50, alt_depth=20,
                filter_status="PASS", expression_tpm=85.0 - (i * 15),
            ),
            peptide=MutantPeptide(
                gene=genes[i], variant_id=f"chr1:{100+i}:G>A",
                mutation=f"p.V{i}I", wt_sequence="A" * 9,
                mut_sequence=peptides[i], position=2,
                length=9, mhc_class=1,
            ),
            binding=BindingPrediction(
                peptide_sequence=peptides[i], allele="DLA-88*001:01",
                affinity_nm=affinities[i], percentile_rank=float(i + 1),
                tool="iedb", mhc_class=1,
            ),
            expression_tpm=85.0 - (i * 15),
            composite_score=0.9 - (i * 0.1),
            rank=i + 1,
            score_components={
                "binding": 0.9 - (i * 0.1),
                "expression": 0.8 - (i * 0.1),
                "vaf": 0.7,
                "self_difference": 0.5,
                "agretopicity": 0.5,
                "caller_agreement": 0.33,
            },
        ))
    return candidates


class TestScoreDistributionChart:

    def test_returns_plotly_figure(self):
        import plotly.graph_objects as go

        from dogneo.ui.charts import score_distribution_chart

        fig = score_distribution_chart(_make_candidates())
        assert isinstance(fig, go.Figure)

    def test_empty_candidates(self):
        import plotly.graph_objects as go

        from dogneo.ui.charts import score_distribution_chart

        fig = score_distribution_chart([])
        assert isinstance(fig, go.Figure)


class TestBindingHeatmap:

    def test_returns_plotly_figure(self):
        import plotly.graph_objects as go

        from dogneo.ui.charts import binding_heatmap

        fig = binding_heatmap(_make_candidates())
        assert isinstance(fig, go.Figure)


class TestRadarChart:

    def test_returns_plotly_figure(self):
        import plotly.graph_objects as go

        from dogneo.ui.charts import score_radar_chart

        candidate = _make_candidates(1)[0]
        fig = score_radar_chart(candidate)
        assert isinstance(fig, go.Figure)


class TestScatterPlot:

    def test_returns_plotly_figure(self):
        import plotly.graph_objects as go

        from dogneo.ui.charts import candidate_scatter

        fig = candidate_scatter(_make_candidates(), x="binding_affinity_nm", y="expression_tpm")
        assert isinstance(fig, go.Figure)
