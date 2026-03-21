"""Multi-factor neoantigen candidate scoring and ranking.

Combines binding affinity, expression level, variant allele frequency,
self-similarity, and other factors into a composite immunogenicity score.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.variants import SomaticVariant

logger = logging.getLogger(__name__)


@dataclass
class NeoantigenCandidate:
    """A scored neoantigen candidate combining all analysis layers.

    Attributes:
        variant: Source somatic variant.
        peptide: Mutant peptide.
        binding: Best binding prediction for this peptide.
        expression_tpm: Gene expression level (TPM).
        composite_score: Final immunogenicity score (higher = better).
        rank: Final rank position.
        score_components: Individual score components for transparency.
    """

    variant: SomaticVariant
    peptide: MutantPeptide
    binding: BindingPrediction
    expression_tpm: float = 0.0
    composite_score: float = 0.0
    rank: int = 0
    score_components: dict[str, float] = field(default_factory=dict)

    @property
    def candidate_id(self) -> str:
        """Unique identifier for this candidate."""
        return f"{self.variant.variant_id}|{self.peptide.mut_sequence}|{self.binding.allele}"

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "rank": self.rank,
            "gene": self.variant.gene,
            "variant_id": self.variant.variant_id,
            "mutation": self.peptide.mutation,
            "mutant_peptide": self.peptide.mut_sequence,
            "wildtype_peptide": self.peptide.wt_sequence,
            "peptide_length": self.peptide.length,
            "mhc_class": self.peptide.mhc_class,
            "allele": self.binding.allele,
            "binding_affinity_nm": self.binding.affinity_nm,
            "binding_percentile": self.binding.percentile_rank,
            "binding_tool": self.binding.tool,
            "vaf": self.variant.vaf,
            "expression_tpm": self.expression_tpm,
            "composite_score": round(self.composite_score, 4),
            "score_components": {k: round(v, 4) for k, v in self.score_components.items()},
            "num_callers": len(self.variant.callers),
        }


# ---------------------------------------------------------------------------
# Scoring Functions
# ---------------------------------------------------------------------------

@dataclass
class ScoringWeights:
    """Weights for multi-factor scoring.

    Defaults represent a balanced scoring approach. Adjust based on
    experimental validation or domain expertise.
    """

    binding_affinity: float = 0.30
    expression: float = 0.20
    vaf: float = 0.15
    self_difference: float = 0.15
    agretopicity: float = 0.10
    caller_agreement: float = 0.10


def _score_binding(affinity_nm: float) -> float:
    """Score binding affinity (0-1, higher = stronger binding).

    Uses log-linear transform: 50 nM → ~1.0, 500 nM → ~0.5, 5000 nM → ~0.0
    """
    if affinity_nm <= 0:
        return 1.0
    score = 1.0 - (math.log10(affinity_nm) - math.log10(50)) / (math.log10(5000) - math.log10(50))
    return max(0.0, min(1.0, score))


def _score_expression(tpm: float) -> float:
    """Score expression level (0-1, higher = more expressed).

    Uses log-linear transform: TPM 1 → ~0.0, TPM 100 → ~1.0
    """
    if tpm <= 0:
        return 0.0
    score = math.log10(tpm + 1) / math.log10(101)
    return max(0.0, min(1.0, score))


def _score_vaf(vaf: float) -> float:
    """Score variant allele frequency (0-1, higher = more clonal).

    Linear: VAF 0 → 0, VAF 0.5+ → 1.0
    """
    return min(1.0, vaf * 2.0)


def _score_self_difference(wt_sequence: str, mut_sequence: str) -> float:
    """Score difference between wild-type and mutant peptide (0-1).

    More differences = higher score (more immunogenic potential).
    """
    if len(wt_sequence) != len(mut_sequence):
        return 1.0  # Length change = high difference

    differences = sum(1 for a, b in zip(wt_sequence, mut_sequence) if a != b)
    return min(1.0, differences / 3.0)  # Cap at 3 differences


def _score_agretopicity(wt_binding_nm: float, mut_binding_nm: float) -> float:
    """Score agretopicity (ratio of wt/mut binding).

    Higher agretopicity = mutant binds much better than wild-type.
    Agretopicity = wt_affinity / mut_affinity. Score > 1 means mut binds better.
    """
    if mut_binding_nm <= 0:
        return 1.0
    if wt_binding_nm <= 0:
        wt_binding_nm = 50000  # Assume very weak WT binding

    ratio = wt_binding_nm / mut_binding_nm
    # Normalize: ratio 1 → 0.5, ratio 10 → ~1.0, ratio 0.1 → ~0.0
    score = 0.5 + 0.5 * math.tanh(math.log10(ratio))
    return max(0.0, min(1.0, score))


def _score_caller_agreement(num_callers: int, max_callers: int = 3) -> float:
    """Score based on number of variant callers that agree."""
    return min(1.0, num_callers / max_callers)


# ---------------------------------------------------------------------------
# Main Ranking Function
# ---------------------------------------------------------------------------

def rank_candidates(
    candidates: list[NeoantigenCandidate],
    weights: ScoringWeights | None = None,
) -> list[NeoantigenCandidate]:
    """Score and rank neoantigen candidates by composite immunogenicity.

    Args:
        candidates: List of unscored candidates.
        weights: Scoring weights. Uses defaults if None.

    Returns:
        Same list sorted by composite_score (descending), with scores and ranks set.
    """
    if weights is None:
        weights = ScoringWeights()

    for candidate in candidates:
        components: dict[str, float] = {}

        # Binding affinity score
        components["binding"] = _score_binding(candidate.binding.affinity_nm)

        # Expression score
        components["expression"] = _score_expression(candidate.expression_tpm)

        # VAF score
        components["vaf"] = _score_vaf(candidate.variant.vaf)

        # Self-difference score
        components["self_difference"] = _score_self_difference(
            candidate.peptide.wt_sequence,
            candidate.peptide.mut_sequence,
        )

        # Agretopicity (placeholder: would need WT binding prediction)
        components["agretopicity"] = 0.5  # Default neutral

        # Caller agreement
        components["caller_agreement"] = _score_caller_agreement(
            len(candidate.variant.callers)
        )

        # Composite score
        composite = (
            weights.binding_affinity * components["binding"]
            + weights.expression * components["expression"]
            + weights.vaf * components["vaf"]
            + weights.self_difference * components["self_difference"]
            + weights.agretopicity * components["agretopicity"]
            + weights.caller_agreement * components["caller_agreement"]
        )

        candidate.composite_score = composite
        candidate.score_components = components

    # Sort descending by score
    candidates.sort(key=lambda c: c.composite_score, reverse=True)

    # Assign ranks
    for i, candidate in enumerate(candidates, start=1):
        candidate.rank = i

    logger.info("Ranked %d neoantigen candidates", len(candidates))
    return candidates


def build_candidates(
    variants: list[SomaticVariant],
    peptides_by_variant: dict[str, list[MutantPeptide]],
    predictions_by_peptide: dict[str, list[BindingPrediction]],
) -> list[NeoantigenCandidate]:
    """Assemble NeoantigenCandidate objects from computed components.

    Args:
        variants: Somatic variants with expression annotated.
        peptides_by_variant: Peptides keyed by variant_id.
        predictions_by_peptide: Binding predictions keyed by peptide sequence.

    Returns:
        List of NeoantigenCandidate objects (unscored).
    """
    candidates: list[NeoantigenCandidate] = []

    for variant in variants:
        peptides = peptides_by_variant.get(variant.variant_id, [])
        for peptide in peptides:
            preds = predictions_by_peptide.get(peptide.mut_sequence, [])
            for pred in preds:
                candidates.append(NeoantigenCandidate(
                    variant=variant,
                    peptide=peptide,
                    binding=pred,
                    expression_tpm=variant.expression_tpm,
                ))

    logger.info("Built %d neoantigen candidates", len(candidates))
    return candidates
