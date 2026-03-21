"""Rank pipeline service layer.

Extracts the core ranking logic from cli.py into a reusable service
that can be called by CLI, UI, notebooks, or tests.
"""
from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import ProteinDatabase, generate_peptides
from dogneo.core.ranking import NeoantigenCandidate, build_candidates, rank_candidates
from dogneo.core.variants import load_vcf, filter_variants
from dogneo.export.exporters import export_fasta, export_json, export_tsv

logger = logging.getLogger(__name__)


@dataclass
class RankInput:
    """Input parameters for the rank pipeline."""

    vcf_path: Path
    expression_path: Path | None = None
    sample_id: str = "SAMPLE"
    alleles: list[str] = field(default_factory=list)
    mhci_lengths: list[int] = field(default_factory=lambda: [8, 9, 10, 11])
    protein_db_path: Path | None = None
    binding_tool: str = "auto"
    llm_tier: str = "none"
    formats: list[str] = field(default_factory=lambda: ["tsv", "json"])


@dataclass
class RankResult:
    """Output of the rank pipeline."""

    candidates: list[NeoantigenCandidate]
    variants_total: int
    variants_coding: int
    peptides_total: int
    binding_tool_used: str
    explanations: dict[str, str] = field(default_factory=dict)
    alleles_used: list[str] = field(default_factory=list)


def _auto_load_alleles() -> list[str]:
    """Load bundled DLA alleles if available."""
    from dogneo.data.manager import ReferenceDataManager

    mgr = ReferenceDataManager()
    dla_path = mgr.get_dla_alleles_path()
    if dla_path.exists():
        alleles = [
            line.strip()
            for line in open(dla_path)
            if line.strip() and not line.strip().startswith("#")
        ]
        return alleles[:4]
    return []


def _resolve_binding_tool(tool: str) -> str:
    """Resolve binding tool with auto-fallback.

    Priority: netmhcpan (local) → iedb (remote) → none
    """
    if tool != "auto":
        return tool
    if shutil.which("netMHCpan"):
        return "netmhcpan"
    return "iedb"


def _auto_detect_proteome(explicit_path: Path | None) -> Path | None:
    """Resolve proteome path: explicit > cached > None."""
    if explicit_path:
        return explicit_path

    from dogneo.data.manager import ReferenceDataManager

    mgr = ReferenceDataManager()
    cached = mgr.get_proteome_path()
    return cached


def run_rank_pipeline(inp: RankInput, output_dir: Path) -> RankResult:
    """Execute the full rank pipeline.

    Args:
        inp: Pipeline input parameters.
        output_dir: Base output directory.

    Returns:
        RankResult with candidates and metadata.
    """
    sample_dir = output_dir / inp.sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Auto-load alleles if none provided
    alleles = inp.alleles
    if not alleles:
        alleles = _auto_load_alleles()
        if alleles:
            logger.info("Auto-loaded %d DLA alleles from bundled data", len(alleles))

    # Step 1: Load and filter variants
    variants = load_vcf(str(inp.vcf_path))
    coding = filter_variants(variants)

    # Step 2: Load expression (optional)
    if inp.expression_path:
        from dogneo.core.expression import annotate_expression, load_expression

        expr_data = load_expression(str(inp.expression_path))
        coding = annotate_expression(coding, expr_data)

    # Step 3: Load proteome and generate peptides
    pdb = ProteinDatabase()
    proteome_path = _auto_detect_proteome(inp.protein_db_path)
    if proteome_path:
        pdb.load_fasta(str(proteome_path))

    peptides_by_variant: dict[str, list] = {}
    for v in coding:
        peps = generate_peptides(v, pdb, lengths=inp.mhci_lengths)
        if peps:
            peptides_by_variant[v.variant_id] = peps

    total_peptides = sum(len(peps) for peps in peptides_by_variant.values())

    # Step 4: Binding prediction + ranking
    binding_tool_used = _resolve_binding_tool(inp.binding_tool)
    candidates_list: list[NeoantigenCandidate] = []

    # Collect all unique peptide sequences for batch prediction
    all_peptides: list[str] = []
    if peptides_by_variant:
        for peps in peptides_by_variant.values():
            for p in peps:
                if p.mut_sequence not in all_peptides:
                    all_peptides.append(p.mut_sequence)

    # Run binding prediction if tool available
    predictions_by_key: dict[tuple[str, str], BindingPrediction] = {}
    if binding_tool_used == "iedb" and all_peptides and alleles:
        try:
            from dogneo.core.iedb_client import IEDBClient

            client = IEDBClient()
            preds = client.predict_batch(all_peptides, alleles)
            for pred in preds:
                predictions_by_key[(pred.peptide_sequence, pred.allele)] = pred
            logger.info("IEDB returned %d binding predictions", len(preds))
        except Exception as e:
            logger.warning("IEDB binding prediction failed: %s", e)
            binding_tool_used = "none"

    # Build candidates with binding predictions (or unscored)
    if peptides_by_variant:
        for variant in coding:
            peptides = peptides_by_variant.get(variant.variant_id, [])
            for peptide in peptides:
                for allele in alleles or ["DLA-88*unknown"]:
                    key = (peptide.mut_sequence, allele)
                    if key in predictions_by_key:
                        binding = predictions_by_key[key]
                    else:
                        binding = BindingPrediction(
                            peptide_sequence=peptide.mut_sequence,
                            allele=allele,
                            affinity_nm=float("nan"),
                            percentile_rank=float("nan"),
                            tool="pending",
                            mhc_class=1,
                        )
                    candidates_list.append(
                        NeoantigenCandidate(
                            variant=variant,
                            peptide=peptide,
                            binding=binding,
                            expression_tpm=variant.expression_tpm,
                        )
                    )

    # Rank scored candidates; leave unscored as-is
    if predictions_by_key and candidates_list:
        ranked = rank_candidates(candidates_list)
    else:
        ranked = candidates_list

    # Step 5: Export
    if "tsv" in inp.formats:
        export_tsv(ranked, sample_dir / "candidates.tsv")
    if "json" in inp.formats:
        export_json(ranked, sample_dir / "candidates.json", inp.sample_id)
    if "fasta" in inp.formats:
        export_fasta(ranked, sample_dir / "candidates.fasta")

    return RankResult(
        candidates=ranked,
        variants_total=len(variants),
        variants_coding=len(coding),
        peptides_total=total_peptides,
        binding_tool_used=binding_tool_used,
        alleles_used=alleles,
    )
