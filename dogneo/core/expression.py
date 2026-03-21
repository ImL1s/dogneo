"""RNA-seq expression quantification loading and variant annotation.

Parses output from Salmon or Kallisto to attach gene-level expression
values (TPM) to somatic variants.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

from dogneo.core.variants import SomaticVariant

logger = logging.getLogger(__name__)


@dataclass
class ExpressionData:
    """Gene and transcript-level expression values.

    Attributes:
        gene_tpm: Gene-level TPM values.
        transcript_tpm: Transcript-level TPM values.
        tool: Quantification tool used (salmon, kallisto).
    """

    gene_tpm: dict[str, float] = field(default_factory=dict)
    transcript_tpm: dict[str, float] = field(default_factory=dict)
    tool: str = "unknown"

    def get_gene_tpm(self, gene: str) -> float:
        """Get TPM for a gene, returning 0 if not found."""
        return self.gene_tpm.get(gene, 0.0)

    def get_transcript_tpm(self, transcript_id: str) -> float:
        """Get TPM for a transcript, returning 0 if not found."""
        return self.transcript_tpm.get(transcript_id, 0.0)


def load_salmon_quant(path: str | Path) -> ExpressionData:
    """Load expression data from Salmon quant.sf output.

    Salmon quant.sf columns: Name, Length, EffectiveLength, TPM, NumReads

    Args:
        path: Path to Salmon quant.sf file.

    Returns:
        ExpressionData with transcript-level TPM.
    """
    path = Path(path)
    data = ExpressionData(tool="salmon")

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            transcript_id = row["Name"]
            tpm = float(row["TPM"])
            data.transcript_tpm[transcript_id] = tpm

    logger.info("Loaded %d transcript TPMs from Salmon: %s",
                len(data.transcript_tpm), path.name)
    return data


def load_kallisto_abundance(path: str | Path) -> ExpressionData:
    """Load expression data from Kallisto abundance.tsv output.

    Kallisto columns: target_id, length, eff_length, est_counts, tpm

    Args:
        path: Path to Kallisto abundance.tsv file.

    Returns:
        ExpressionData with transcript-level TPM.
    """
    path = Path(path)
    data = ExpressionData(tool="kallisto")

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            transcript_id = row["target_id"]
            tpm = float(row["tpm"])
            data.transcript_tpm[transcript_id] = tpm

    logger.info("Loaded %d transcript TPMs from Kallisto: %s",
                len(data.transcript_tpm), path.name)
    return data


def load_expression(
    path: str | Path,
    tool: str = "salmon",
) -> ExpressionData:
    """Load expression data from Salmon or Kallisto output.

    Args:
        path: Path to quantification output file.
        tool: Quantification tool name ("salmon" or "kallisto").

    Returns:
        ExpressionData instance.

    Raises:
        ValueError: If tool is not recognized.
    """
    if tool == "salmon":
        return load_salmon_quant(path)
    elif tool == "kallisto":
        return load_kallisto_abundance(path)
    else:
        raise ValueError(f"Unknown expression tool: {tool}. Use 'salmon' or 'kallisto'.")


def aggregate_to_gene_level(
    expression: ExpressionData,
    tx2gene: dict[str, str] | None = None,
) -> ExpressionData:
    """Aggregate transcript-level TPM to gene-level.

    Args:
        expression: Transcript-level expression data.
        tx2gene: Optional transcript-to-gene mapping. If None, uses
            transcript IDs as gene names (truncating version suffixes).

    Returns:
        Same ExpressionData with gene_tpm populated.
    """
    gene_sums: dict[str, float] = {}

    for tx_id, tpm in expression.transcript_tpm.items():
        if tx2gene is not None:
            gene = tx2gene.get(tx_id, tx_id)
        else:
            # Strip version suffix (e.g., ENSCAFT00000012345.3 → ENSCAFT00000012345)
            gene = tx_id.split(".")[0]

        gene_sums[gene] = gene_sums.get(gene, 0.0) + tpm

    expression.gene_tpm = gene_sums
    logger.info("Aggregated to %d gene-level TPMs", len(gene_sums))
    return expression


def annotate_expression(
    variants: list[SomaticVariant],
    expression: ExpressionData,
) -> list[SomaticVariant]:
    """Attach expression values (TPM) to each variant's gene.

    Args:
        variants: List of somatic variants.
        expression: Expression data with gene-level TPM.

    Returns:
        Same variant list with expression_tpm populated.
    """
    annotated = 0
    for v in variants:
        if v.gene:
            tpm = expression.get_gene_tpm(v.gene)
            if tpm == 0.0 and v.transcript_id:
                tpm = expression.get_transcript_tpm(v.transcript_id)
            v.expression_tpm = tpm
            if tpm > 0:
                annotated += 1

    logger.info("Annotated expression for %d / %d variants", annotated, len(variants))
    return variants
