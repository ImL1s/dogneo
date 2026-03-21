"""Plotly chart builders for DogNeo UI.

All charts use plotly_dark template for consistent dark-theme styling.
"""
from __future__ import annotations

import plotly.graph_objects as go

from dogneo.core.ranking import NeoantigenCandidate

TEMPLATE = "plotly_dark"


def score_distribution_chart(candidates: list[NeoantigenCandidate]) -> go.Figure:
    """Bar chart of composite scores by rank.

    Args:
        candidates: Ranked neoantigen candidates.

    Returns:
        Plotly Figure.
    """
    if not candidates:
        fig = go.Figure()
        fig.update_layout(
            template=TEMPLATE,
            title="Score Distribution (no candidates)",
        )
        return fig

    genes = [f"{c.variant.gene} {c.peptide.mutation}" for c in candidates]
    scores = [c.composite_score for c in candidates]
    affinities = [c.binding.affinity_nm for c in candidates]

    # Color by binding strength
    colors = []
    for a in affinities:
        if a < 50:
            colors.append("#4ade80")  # green — strong
        elif a < 500:
            colors.append("#facc15")  # yellow — weak
        else:
            colors.append("#94a3b8")  # grey — non-binder

    fig = go.Figure(go.Bar(
        x=list(range(1, len(candidates) + 1)),
        y=scores,
        marker_color=colors,
        text=genes,
        hovertemplate=(
            "<b>%{text}</b><br>"
            "Rank: %{x}<br>"
            "Score: %{y:.4f}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        template=TEMPLATE,
        title="Neoantigen Candidate Scores",
        xaxis_title="Rank",
        yaxis_title="Composite Score",
        showlegend=False,
    )
    return fig


def binding_heatmap(candidates: list[NeoantigenCandidate]) -> go.Figure:
    """Heatmap of binding affinity: peptide × allele.

    Args:
        candidates: Ranked neoantigen candidates.

    Returns:
        Plotly Figure.
    """
    if not candidates:
        fig = go.Figure()
        fig.update_layout(template=TEMPLATE, title="Binding Heatmap (no data)")
        return fig

    # Collect unique peptides and alleles
    peptides: list[str] = []
    alleles: list[str] = []
    affinity_map: dict[tuple[str, str], float] = {}

    for c in candidates:
        pep_label = f"{c.variant.gene} {c.peptide.mut_sequence}"
        if pep_label not in peptides:
            peptides.append(pep_label)
        if c.binding.allele not in alleles:
            alleles.append(c.binding.allele)
        affinity_map[(pep_label, c.binding.allele)] = c.binding.affinity_nm

    # Build z-matrix
    z = []
    for pep in peptides:
        row = [affinity_map.get((pep, allele), float("nan")) for allele in alleles]
        z.append(row)

    fig = go.Figure(go.Heatmap(
        z=z,
        x=alleles,
        y=peptides,
        colorscale="RdYlGn_r",  # Red=strong binding (low nM), Green=weak
        colorbar_title="Affinity (nM)",
        hovertemplate=(
            "Peptide: %{y}<br>"
            "Allele: %{x}<br>"
            "Affinity: %{z:.1f} nM<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        template=TEMPLATE,
        title="Peptide–DLA Binding Affinity",
        xaxis_title="DLA Allele",
        yaxis_title="Peptide",
    )
    return fig


def score_radar_chart(candidate: NeoantigenCandidate) -> go.Figure:
    """Radar chart of 6 scoring components for a single candidate.

    Args:
        candidate: A single scored candidate.

    Returns:
        Plotly Figure.
    """
    components = candidate.score_components
    categories = list(components.keys())
    values = list(components.values())

    # Close the polygon
    categories.append(categories[0])
    values.append(values[0])

    fig = go.Figure(go.Scatterpolar(
        r=values,
        theta=categories,
        fill="toself",
        fillcolor="rgba(56, 189, 248, 0.3)",
        line_color="#38bdf8",
        name=f"{candidate.variant.gene} {candidate.peptide.mutation}",
    ))
    fig.update_layout(
        template=TEMPLATE,
        title=f"Score Components — {candidate.variant.gene} {candidate.peptide.mutation}",
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 1]),
        ),
        showlegend=False,
    )
    return fig


def candidate_scatter(
    candidates: list[NeoantigenCandidate],
    x: str = "binding_affinity_nm",
    y: str = "expression_tpm",
) -> go.Figure:
    """Configurable scatter plot of candidate properties.

    Args:
        candidates: Ranked candidates.
        x: X-axis field from candidate.to_dict().
        y: Y-axis field from candidate.to_dict().

    Returns:
        Plotly Figure.
    """
    if not candidates:
        fig = go.Figure()
        fig.update_layout(template=TEMPLATE, title="Scatter (no data)")
        return fig

    dicts = [c.to_dict() for c in candidates]
    x_vals = [d.get(x, 0) for d in dicts]
    y_vals = [d.get(y, 0) for d in dicts]
    labels = [f"{d['gene']} {d['mutation']}" for d in dicts]
    scores = [d["composite_score"] for d in dicts]

    fig = go.Figure(go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="markers",
        marker=dict(
            size=10,
            color=scores,
            colorscale="Viridis",
            colorbar_title="Score",
        ),
        text=labels,
        hovertemplate=(
            "<b>%{text}</b><br>"
            f"{x}: %{{x:.1f}}<br>"
            f"{y}: %{{y:.1f}}<br>"
            "Score: %{marker.color:.4f}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        template=TEMPLATE,
        title=f"{y} vs {x}",
        xaxis_title=x.replace("_", " ").title(),
        yaxis_title=y.replace("_", " ").title(),
    )
    return fig
