"""Pipeline step-by-step explainer using LLM.

Generates plain-language explanations at each pipeline stage,
making results accessible to both researchers and dog owners.
Falls back to structured summaries when no LLM is available.
"""
from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

from dogneo.core.variants import SomaticVariant
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.llm.router import LLMRouter, TaskType

if TYPE_CHECKING:
    from dogneo.app.rank_pipeline import RankResult

logger = logging.getLogger(__name__)


class PipelineExplainer:
    """Generate plain-language explanations at each pipeline step.

    When an LLM router is available, uses it for rich narrative.
    When not available, returns structured text summaries.

    Args:
        router: Optional LLM router for AI-powered explanations.
    """

    def __init__(self, router: LLMRouter | None = None):
        self.router = router

    def _generate(self, prompt: str, fallback: str) -> str:
        """Try LLM generation, fall back to static text on failure."""
        if self.router is None:
            return fallback
        try:
            return self.router.generate(prompt, task_type=TaskType.SUMMARIZE)
        except Exception as e:
            logger.warning("LLM explanation failed: %s", e)
            return fallback

    def explain_variants(
        self, variants: list[SomaticVariant], coding: list[SomaticVariant]
    ) -> str:
        """Explain variant loading and filtering results."""
        genes = sorted(set(v.gene for v in coding if v.gene))
        gene_str = ", ".join(genes) if genes else "unknown genes"

        fallback = (
            f"Found {len(variants)} total variants, {len(coding)} affect proteins. "
            f"Genes with coding mutations: {gene_str}."
        )

        prompt = (
            f"You are a veterinary genomics assistant explaining canine tumor analysis results. "
            f"Explain in plain language for a dog owner:\n\n"
            f"We analyzed the tumor DNA and found {len(variants)} total mutations. "
            f"Of these, {len(coding)} affect protein sequences (coding mutations). "
            f"The affected genes are: {gene_str}.\n\n"
            f"Briefly explain what this means for the dog's cancer analysis. "
            f"Keep it under 100 words, compassionate tone, no jargon."
        )

        return self._generate(prompt, fallback)

    def explain_peptides(
        self, peptides_by_variant: dict[str, list], total: int
    ) -> str:
        """Explain peptide generation results."""
        n_variants = len(peptides_by_variant)
        fallback = (
            f"Generated {total} candidate peptides from {n_variants} mutations. "
            f"Each peptide is a short protein fragment (8-11 amino acids) "
            f"that the immune system might recognize."
        )

        prompt = (
            f"You are a veterinary genomics assistant. Explain to a dog owner:\n\n"
            f"From {n_variants} mutations, we generated {total} candidate peptide fragments. "
            f"These are short pieces of mutated protein (8-11 amino acids long). "
            f"The dog's immune system needs to recognize these fragments to fight the tumor.\n\n"
            f"Explain what happens next (binding prediction). Under 80 words, simple language."
        )

        return self._generate(prompt, fallback)

    def explain_binding(self, candidates: list[NeoantigenCandidate]) -> str:
        """Explain binding prediction results."""
        strong = sum(1 for c in candidates if not math.isnan(c.binding.affinity_nm) and c.binding.affinity_nm < 50)
        weak = sum(1 for c in candidates if not math.isnan(c.binding.affinity_nm) and 50 <= c.binding.affinity_nm < 500)
        total = len(candidates)

        fallback = (
            f"Evaluated {total} peptide-allele combinations. "
            f"{strong} are strong binders (<50 nM), {weak} are weak binders (<500 nM). "
            f"Strong binders are the best vaccine candidates."
        )

        prompt = (
            f"You are a veterinary genomics assistant. Explain to a dog owner:\n\n"
            f"We tested {total} peptide fragments against your dog's immune system markers (DLA alleles). "
            f"{strong} bind very strongly (excellent vaccine candidates), "
            f"{weak} bind moderately (possible candidates).\n\n"
            f"Explain what binding strength means for vaccine effectiveness. "
            f"Under 80 words, compassionate, no jargon."
        )

        return self._generate(prompt, fallback)

    def explain_ranking(self, top_candidates: list[NeoantigenCandidate]) -> str:
        """Explain the top ranked candidates."""
        if not top_candidates:
            return "No candidates were ranked."

        top = top_candidates[0]
        details = []
        for i, c in enumerate(top_candidates[:3], 1):
            details.append(
                f"#{i}: {c.variant.gene} {c.peptide.mutation} "
                f"(binding: {c.binding.affinity_nm:.0f} nM, "
                f"expression: {c.expression_tpm:.1f} TPM, "
                f"score: {c.composite_score:.3f})"
            )
        details_str = "\n".join(details)

        fallback = (
            f"Top candidate: {top.variant.gene} {top.peptide.mutation} "
            f"with composite score {top.composite_score:.3f}.\n{details_str}"
        )

        prompt = (
            f"You are a veterinary genomics assistant. Explain to a dog owner:\n\n"
            f"Top vaccine candidates ranked by immunogenicity score:\n{details_str}\n\n"
            f"Explain why #{1} ranks highest and what the scores mean. "
            f"Mention the gene's role in cancer if known. "
            f"Under 100 words, compassionate tone."
        )

        return self._generate(prompt, fallback)

    def explain_for_owner(self, result: "RankResult") -> str:
        """Generate a complete narrative summary for dog owners."""
        top_genes = sorted(set(
            c.variant.gene for c in result.candidates[:5] if c.variant.gene
        ))

        fallback = (
            f"Analysis complete. We examined {result.variants_total} mutations in your dog's tumor DNA. "
            f"{result.variants_coding} affect proteins, generating {result.peptides_total} candidate peptides. "
            f"Binding predictions used: {result.binding_tool_used}. "
            f"{len(result.candidates)} candidates were evaluated. "
            f"Top genes: {', '.join(top_genes) if top_genes else 'N/A'}. "
            f"These results should be reviewed by a veterinary oncologist."
        )

        prompt = (
            f"You are a veterinary genomics assistant writing a summary for a dog owner "
            f"whose pet has cancer. Write a compassionate, clear summary:\n\n"
            f"- Total mutations found: {result.variants_total}\n"
            f"- Protein-affecting mutations: {result.variants_coding}\n"
            f"- Candidate peptides generated: {result.peptides_total}\n"
            f"- Binding tool used: {result.binding_tool_used}\n"
            f"- Total candidates evaluated: {len(result.candidates)}\n"
            f"- Top genes involved: {', '.join(top_genes)}\n\n"
            f"Explain what was found, what it means, and recommended next steps "
            f"(consult veterinary oncologist, consider mRNA vaccine synthesis). "
            f"Under 150 words. End with the RUO disclaimer."
        )

        return self._generate(prompt, fallback)

    def literature_check(self, genes: list[str]) -> str:
        """Cross-reference genes against known canine cancer literature."""
        known_canine_genes = {
            "TP53": "Tumor suppressor, frequently mutated in canine osteosarcoma and mammary tumors",
            "BRAF": "Growth signal kinase, V595E mutation common in canine urothelial carcinoma",
            "KRAS": "Oncogene, mutations drive uncontrolled cell growth",
            "PIK3CA": "PI3K pathway, involved in cell survival and proliferation",
            "PTEN": "Tumor suppressor, loss associated with aggressive canine tumors",
            "KIT": "Receptor tyrosine kinase, mutations drive canine mast cell tumors",
        }

        found = {g: known_canine_genes[g] for g in genes if g in known_canine_genes}
        unknown = [g for g in genes if g not in known_canine_genes]

        lines = []
        for gene, desc in found.items():
            lines.append(f"  {gene}: {desc}")

        fallback = "Known canine cancer genes:\n" + "\n".join(lines) if lines else ""
        if unknown:
            fallback += f"\nGenes with limited canine cancer data: {', '.join(unknown)}"

        if not self.router:
            return fallback

        prompt = (
            f"You are a veterinary oncology researcher. For these genes found in a dog's tumor, "
            f"provide a brief literature summary of their known roles in canine cancer:\n\n"
            f"Genes: {', '.join(genes)}\n\n"
            f"For each gene, mention: known canine cancer associations, mutation frequency, "
            f"and any relevant therapeutic implications. Under 150 words."
        )

        return self._generate(prompt, fallback)
