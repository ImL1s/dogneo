"""Export ranked neoantigen candidates to various formats.

Supported formats:
- FASTA: standard peptide sequences for wet lab / synthesis
- TSV: tabular data for spreadsheet analysis
- JSON: structured data with full metadata
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

from dogneo import __version__, RUO_DISCLAIMER
from dogneo.core.ranking import NeoantigenCandidate

logger = logging.getLogger(__name__)


def export_fasta(
    candidates: list[NeoantigenCandidate],
    output_path: str | Path,
    top_n: int | None = None,
) -> Path:
    """Export candidate peptides in FASTA format.

    FASTA headers include gene, mutation, allele, rank, and score.
    Useful for peptide synthesis orders and downstream wet lab work.

    Args:
        candidates: Ranked candidates to export.
        output_path: Output file path.
        top_n: Limit to top N candidates (None for all).

    Returns:
        Path to written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    subset = candidates[:top_n] if top_n else candidates

    with open(output_path, "w", encoding="utf-8") as f:
        # Header comment
        f.write(f"; DogNeo v{__version__} — {RUO_DISCLAIMER}\n")
        f.write(f"; Generated: {datetime.now(timezone.utc).isoformat()}\n")
        f.write(f"; Candidates: {len(subset)}\n\n")

        for c in subset:
            d = c.to_dict()
            header = (
                f">{d['gene']}_{d['mutation']}"
                f" rank={d['rank']}"
                f" allele={d['allele']}"
                f" affinity_nm={d['binding_affinity_nm']:.1f}"
                f" score={d['composite_score']:.4f}"
                f" tpm={d['expression_tpm']:.1f}"
            )
            seq = d["mutant_peptide"]
            f.write(f"{header}\n{seq}\n")

    logger.info("FASTA exported: %s (%d candidates)", output_path, len(subset))
    return output_path


def export_tsv(
    candidates: list[NeoantigenCandidate],
    output_path: str | Path,
    top_n: int | None = None,
) -> Path:
    """Export candidates as tab-separated values.

    Columns: rank, gene, mutation, mutant_peptide, wildtype_peptide,
    allele, binding_affinity_nm, percentile_rank, expression_tpm,
    vaf, composite_score, mhc_class, variant_id.

    Args:
        candidates: Ranked candidates.
        output_path: Output file path.
        top_n: Limit to top N candidates.

    Returns:
        Path to written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    subset = candidates[:top_n] if top_n else candidates

    fieldnames = [
        "rank", "gene", "mutation", "mutant_peptide", "wildtype_peptide",
        "allele", "binding_affinity_nm", "percentile_rank",
        "expression_tpm", "vaf", "composite_score", "mhc_class", "variant_id",
    ]

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        # Comment header
        f.write(f"# DogNeo v{__version__} — {RUO_DISCLAIMER}\n")

        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for c in subset:
            writer.writerow(c.to_dict())

    logger.info("TSV exported: %s (%d candidates)", output_path, len(subset))
    return output_path


def export_json(
    candidates: list[NeoantigenCandidate],
    output_path: str | Path,
    sample_id: str = "",
    top_n: int | None = None,
) -> Path:
    """Export candidates as structured JSON with metadata.

    JSON structure includes pipeline version, timestamp,
    RUO disclaimer, sample info, and full candidate data.

    Args:
        candidates: Ranked candidates.
        output_path: Output file path.
        sample_id: Sample identifier.
        top_n: Limit to top N candidates.

    Returns:
        Path to written file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    subset = candidates[:top_n] if top_n else candidates

    data = {
        "metadata": {
            "tool": "DogNeo",
            "version": __version__,
            "disclaimer": RUO_DISCLAIMER,
            "generated": datetime.now(timezone.utc).isoformat(),
            "sample_id": sample_id,
            "total_candidates": len(candidates),
            "exported_candidates": len(subset),
        },
        "candidates": [c.to_dict() for c in subset],
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    logger.info("JSON exported: %s (%d candidates)", output_path, len(subset))
    return output_path
