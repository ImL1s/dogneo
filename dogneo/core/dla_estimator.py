"""Built-in DLA binding affinity estimator.

Provides approximate binding predictions for canine DLA-88 alleles
without requiring external tools. Based on published binding motifs
for DLA-88 Class I molecules.

This estimator uses a simplified position-specific scoring approach
derived from known DLA-88 anchor preferences. Results are labeled
as "estimated" and should be validated with NetMHCpan when available.

References:
- Venkataraman et al. (2017) PLOS ONE: DLA-88*501:01 binding motif
- Buckley et al. (2018) J. Immunol.: DLA-88*508:01 peptide binding
"""
from __future__ import annotations

import logging
import math

from dogneo.core.binding import BindingPrediction

logger = logging.getLogger(__name__)

# Amino acid hydrophobicity index (Kyte-Doolittle, normalized 0-1)
_HYDROPHOBICITY: dict[str, float] = {
    "I": 1.00, "V": 0.93, "L": 0.84, "F": 0.61, "C": 0.56,
    "M": 0.42, "A": 0.40, "G": 0.16, "T": 0.12, "S": 0.07,
    "W": 0.07, "Y": 0.04, "P": 0.05, "H": 0.00, "D": 0.00,
    "E": 0.00, "N": 0.00, "Q": 0.00, "K": 0.00, "R": 0.00,
}

# DLA-88 anchor preferences (simplified from published motifs)
# MHC-I typically uses P2 (position 2) and Pω (C-terminal) as primary anchors
# DLA-88 generally prefers hydrophobic residues at these positions
_ANCHOR_WEIGHTS = {
    "P2": 2.0,    # Position 2 — primary anchor
    "Pw": 2.5,    # C-terminal — primary anchor (strongest)
    "P1": 0.5,    # N-terminal — secondary
    "P3": 0.3,    # Position 3 — minor
    "mid": 0.2,   # Middle positions — minor
}

# Preferred residues at key positions for DLA-88
_PREFERRED_P2 = set("LIVMFAY")      # Hydrophobic at P2
_PREFERRED_PW = set("LIVMFAYW")     # Hydrophobic at C-terminal
_PREFERRED_P1 = set("RKHDEY")       # Charged/polar at P1


def _score_position(aa: str, position: str) -> float:
    """Score an amino acid at a given anchor position (0-1, higher=better)."""
    hydro = _HYDROPHOBICITY.get(aa, 0.0)

    if position == "P2":
        return 1.0 if aa in _PREFERRED_P2 else hydro * 0.3
    elif position == "Pw":
        return 1.0 if aa in _PREFERRED_PW else hydro * 0.3
    elif position == "P1":
        return 0.8 if aa in _PREFERRED_P1 else 0.3
    else:
        # Middle positions: slight preference for hydrophobic
        return 0.3 + hydro * 0.4


def _estimate_affinity(peptide: str) -> float:
    """Estimate binding affinity for a peptide to DLA-88.

    Uses a simplified anchor-based scoring model.

    Args:
        peptide: Amino acid sequence (8-11 residues).

    Returns:
        Estimated binding affinity in nM (lower = stronger).
    """
    n = len(peptide)
    if n < 8 or n > 15:
        return 50000.0  # Out of range

    # Score key positions
    score = 0.0
    total_weight = 0.0

    # P1 (N-terminal)
    w = _ANCHOR_WEIGHTS["P1"]
    score += w * _score_position(peptide[0], "P1")
    total_weight += w

    # P2
    w = _ANCHOR_WEIGHTS["P2"]
    score += w * _score_position(peptide[1], "P2")
    total_weight += w

    # P3
    if n >= 9:
        w = _ANCHOR_WEIGHTS["P3"]
        score += w * _score_position(peptide[2], "P3")
        total_weight += w

    # Middle positions
    mid_start = 3 if n >= 9 else 2
    mid_end = n - 1
    for i in range(mid_start, mid_end):
        w = _ANCHOR_WEIGHTS["mid"]
        score += w * _score_position(peptide[i], "mid")
        total_weight += w

    # Pω (C-terminal) — strongest anchor
    w = _ANCHOR_WEIGHTS["Pw"]
    score += w * _score_position(peptide[-1], "Pw")
    total_weight += w

    # Normalize to 0-1
    normalized = score / total_weight if total_weight > 0 else 0.0

    # Convert to nM scale (exponential mapping)
    # normalized 1.0 → ~10 nM (very strong)
    # normalized 0.5 → ~500 nM (borderline)
    # normalized 0.0 → ~50000 nM (non-binder)
    affinity_nm = 10.0 * math.exp((1.0 - normalized) * math.log(5000))

    return round(affinity_nm, 1)


def estimate_binding(
    peptides: list[str],
    alleles: list[str],
) -> list[BindingPrediction]:
    """Estimate DLA binding for peptides using built-in model.

    Args:
        peptides: Peptide sequences.
        alleles: DLA allele names (used for labeling; same model for all DLA-88).

    Returns:
        List of BindingPrediction with tool="dogneo-estimator".
    """
    if not peptides or not alleles:
        return []

    results: list[BindingPrediction] = []

    for peptide in peptides:
        affinity = _estimate_affinity(peptide)
        # Approximate percentile rank from affinity
        # <50 nM → ~0.5%, <500 nM → ~2%, >5000 nM → >10%
        if affinity < 50:
            percentile = affinity / 100.0
        elif affinity < 500:
            percentile = 0.5 + (affinity - 50) / 300.0
        else:
            percentile = 2.0 + math.log10(affinity / 500.0) * 5.0

        for allele in alleles:
            results.append(BindingPrediction(
                peptide_sequence=peptide,
                allele=allele,
                affinity_nm=affinity,
                percentile_rank=round(percentile, 2),
                tool="dogneo-estimator",
                mhc_class=1,
            ))

    logger.info("Estimated binding for %d peptide-allele pairs", len(results))
    return results
