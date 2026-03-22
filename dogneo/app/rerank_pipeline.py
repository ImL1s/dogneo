"""Rerank pipeline — import external binding results and re-score candidates.

Supports importing from NetMHCpan, MHCflurry, or generic TSV formats.
"""
from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from dogneo.core.binding import BindingPrediction
from dogneo.core.ranking import MutantPeptide, NeoantigenCandidate, SomaticVariant, rank_candidates
from dogneo.export.exporters import export_fasta, export_json, export_tsv

logger = logging.getLogger(__name__)


@dataclass
class RerankInput:
    """Input parameters for the rerank pipeline."""

    candidates_path: Path
    binding_path: Path
    binding_format: str = "auto"
    formats: list[str] = field(default_factory=lambda: ["tsv", "json"])


@dataclass
class RerankResult:
    """Output of the rerank pipeline."""

    candidates: list[NeoantigenCandidate]
    binding_matched: int
    binding_unmatched: int


def detect_binding_format(path: Path) -> str:
    """Auto-detect binding result format from file header.

    Returns:
        Format string: 'netmhcpan', 'mhcflurry', or 'tsv'.
    """
    with open(path) as f:
        header_line = f.readline().strip()

    header = header_line.split("\t")

    if "Aff(nM)" in header or ("MHC" in header and "Peptide" in header and "%Rank_BA" in header):
        return "netmhcpan"
    if "mhcflurry_affinity" in header:
        return "mhcflurry"
    return "tsv"


def parse_binding_results(
    path: Path,
    fmt: str = "auto",
) -> dict[tuple[str, str], BindingPrediction]:
    """Parse binding results from external tools.

    Args:
        path: Path to binding results file.
        fmt: Format ('auto', 'netmhcpan', 'mhcflurry', 'tsv').

    Returns:
        Dict mapping (peptide, allele) to BindingPrediction.
    """
    if fmt == "auto":
        fmt = detect_binding_format(path)

    results: dict[tuple[str, str], BindingPrediction] = {}

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                if fmt == "netmhcpan":
                    peptide = row.get("Peptide", "")
                    allele = row.get("MHC", "")
                    affinity = float(row.get("Aff(nM)", row.get("nM", "99999")))
                    percentile = float(row.get("%Rank_BA", row.get("%Rank", "100")))
                    tool = "netmhcpan"
                elif fmt == "mhcflurry":
                    peptide = row.get("peptide", "")
                    allele = row.get("allele", "")
                    affinity = float(row.get("mhcflurry_affinity", "99999"))
                    percentile = float(row.get("mhcflurry_affinity_percentile", "100"))
                    tool = "mhcflurry"
                else:  # generic tsv
                    peptide = row.get("peptide", "")
                    allele = row.get("allele", "")
                    affinity = float(row.get("affinity_nm", "99999"))
                    percentile = float(row.get("percentile_rank", "100"))
                    tool = "external"

                if peptide and allele:
                    results[(peptide, allele)] = BindingPrediction(
                        peptide_sequence=peptide,
                        allele=allele,
                        affinity_nm=affinity,
                        percentile_rank=percentile,
                        tool=tool,
                        mhc_class=1,
                    )
            except (ValueError, KeyError) as e:
                logger.debug("Skipping row: %s", e)
                continue

    logger.info("Parsed %d binding predictions from %s (%s format)", len(results), path.name, fmt)
    return results


def _parse_variant_id(vid: str) -> tuple[str, int, str, str]:
    """Parse variant_id string 'chr:pos:ref>alt' into components."""
    try:
        parts = vid.split(":")
        chrom = parts[0] if len(parts) >= 1 else ""
        pos = int(parts[1]) if len(parts) >= 2 else 0
        allele_part = parts[2] if len(parts) >= 3 else ">"
        ref, alt = allele_part.split(">", 1) if ">" in allele_part else ("", "")
        return chrom, pos, ref, alt
    except (ValueError, IndexError):
        return "", 0, "", ""


def _candidate_from_dict(d: dict) -> NeoantigenCandidate:
    """Reconstruct a NeoantigenCandidate from a serialized dict."""
    chrom, pos, ref, alt = _parse_variant_id(d.get("variant_id", ""))
    variant = SomaticVariant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene=d.get("gene", ""),
        effect="missense_variant",
        hgvs_p=d.get("mutation", ""),
        vaf=float(d.get("vaf", 0)),
        depth=10,
        alt_depth=4,
        filter_status="PASS",
        expression_tpm=float(d.get("expression_tpm", 0)),
    )
    peptide = MutantPeptide(
        gene=d.get("gene", ""),
        variant_id=d.get("variant_id", ""),
        mutation=d.get("mutation", ""),
        wt_sequence=d.get("wildtype_peptide", ""),
        mut_sequence=d.get("mutant_peptide", ""),
        position=0,
        length=int(d.get("peptide_length", 9)),
        mhc_class=int(d.get("mhc_class", 1)),
    )
    binding = BindingPrediction(
        peptide_sequence=d.get("mutant_peptide", ""),
        allele=d.get("allele", ""),
        affinity_nm=float(d.get("binding_affinity_nm", "nan")),
        percentile_rank=float(d.get("binding_percentile", "nan")),
        tool=d.get("binding_tool", "pending"),
        mhc_class=int(d.get("mhc_class", 1)),
    )
    return NeoantigenCandidate(
        variant=variant,
        peptide=peptide,
        binding=binding,
        expression_tpm=float(d.get("expression_tpm", 0)),
    )


def run_rerank_pipeline(inp: RerankInput, output_dir: Path) -> RerankResult:
    """Import external binding results, merge with candidates, and re-score.

    Args:
        inp: Rerank input parameters.
        output_dir: Output directory for re-ranked results.

    Returns:
        RerankResult with re-scored candidates.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load candidates JSON
    with open(inp.candidates_path) as f:
        data = json.load(f)

    candidate_dicts = data.get("candidates", [])
    candidates = [_candidate_from_dict(d) for d in candidate_dicts]

    # Parse binding results
    binding_results = parse_binding_results(inp.binding_path, fmt=inp.binding_format)

    # Merge binding into candidates
    matched = 0
    for candidate in candidates:
        key = (candidate.peptide.mut_sequence, candidate.binding.allele)
        if key in binding_results:
            candidate.binding = binding_results[key]
            matched += 1

    unmatched = len(candidates) - matched
    if unmatched > 0:
        logger.warning("%d candidates had no matching binding prediction", unmatched)

    # Re-score and rank
    ranked = rank_candidates(candidates)

    # Export
    if "tsv" in inp.formats:
        export_tsv(ranked, output_dir / "candidates.tsv")
    if "json" in inp.formats:
        sample_id = data.get("metadata", {}).get("sample_id", "RERANKED")
        export_json(ranked, output_dir / "candidates.json", sample_id)
    if "fasta" in inp.formats:
        export_fasta(ranked, output_dir / "candidates.fasta")

    logger.info("Re-ranked %d candidates (%d with binding data)", len(ranked), matched)

    return RerankResult(
        candidates=ranked,
        binding_matched=matched,
        binding_unmatched=unmatched,
    )
