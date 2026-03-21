"""Somatic variant loading and filtering from annotated VCF files.

Supports VCF files annotated by SnpEff or VEP, extracting protein-level
effects for downstream peptide generation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

logger = logging.getLogger(__name__)


@dataclass
class SomaticVariant:
    """A somatic variant with annotation and quality information.

    Attributes:
        chrom: Chromosome name.
        pos: 1-based genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        gene: Gene symbol (from annotation).
        transcript_id: Transcript identifier.
        effect: Variant effect type (e.g., missense_variant, frameshift).
        hgvs_p: Protein-level HGVS notation (e.g., p.V600E).
        hgvs_c: Coding-level HGVS notation.
        vaf: Variant allele frequency in tumor.
        depth: Total read depth at position.
        alt_depth: Alternate allele read depth.
        filter_status: VCF FILTER field value.
        expression_tpm: Gene expression level (TPM), set during annotation.
        callers: Set of variant callers that detected this variant.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    transcript_id: str = ""
    effect: str = ""
    hgvs_p: str = ""
    hgvs_c: str = ""
    vaf: float = 0.0
    depth: int = 0
    alt_depth: int = 0
    filter_status: str = "PASS"
    expression_tpm: float = 0.0
    callers: set[str] = field(default_factory=set)

    @property
    def variant_id(self) -> str:
        """Unique identifier for this variant."""
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    @property
    def is_coding(self) -> bool:
        """Whether this variant affects protein coding."""
        coding_effects = {
            "missense_variant", "frameshift_variant", "inframe_insertion",
            "inframe_deletion", "stop_gained", "stop_lost",
            "start_lost", "protein_altering_variant",
        }
        return self.effect in coding_effects

    def __hash__(self) -> int:
        return hash(self.variant_id)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SomaticVariant):
            return NotImplemented
        return self.variant_id == other.variant_id


# ---------------------------------------------------------------------------
# VCF Parsing
# ---------------------------------------------------------------------------

def _parse_info_field(info_str: str) -> dict[str, str]:
    """Parse VCF INFO field into key-value dict."""
    info: dict[str, str] = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = "true"
    return info


def _extract_snpeff_annotation(ann_str: str) -> dict[str, str]:
    """Extract relevant fields from SnpEff ANN field.

    SnpEff ANN format (pipe-separated):
    Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|
    Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos|CDS.pos|AA.pos|...
    """
    parts = ann_str.split("|")
    if len(parts) < 11:
        return {}
    return {
        "effect": parts[1],
        "gene": parts[3],
        "transcript_id": parts[6],
        "hgvs_c": parts[9],
        "hgvs_p": parts[10],
    }


def _extract_vep_annotation(csq_str: str) -> dict[str, str]:
    """Extract relevant fields from VEP CSQ field.

    VEP CSQ format varies by configuration. Assumes default:
    Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...
    """
    parts = csq_str.split("|")
    if len(parts) < 7:
        return {}
    return {
        "effect": parts[1],
        "gene": parts[3],
        "transcript_id": parts[6],
        "hgvs_c": parts.get(10, "") if len(parts) > 10 else "",
        "hgvs_p": parts[11] if len(parts) > 11 else "",
    }


def _parse_vcf_line(line: str) -> SomaticVariant | None:
    """Parse a single VCF data line into a SomaticVariant."""
    parts = line.strip().split("\t")
    if len(parts) < 8:
        return None

    chrom, pos_str, _id, ref, alt, _qual, filt, info_str = parts[:8]

    info = _parse_info_field(info_str)

    # Extract annotation (SnpEff ANN or VEP CSQ)
    ann: dict[str, str] = {}
    if "ANN" in info:
        # Take first annotation entry
        first_ann = info["ANN"].split(",")[0]
        ann = _extract_snpeff_annotation(first_ann)
    elif "CSQ" in info:
        first_csq = info["CSQ"].split(",")[0]
        ann = _extract_vep_annotation(first_csq)

    # Extract VAF from tumor sample (FORMAT field)
    vaf = 0.0
    depth = 0
    alt_depth = 0
    if len(parts) > 9:
        fmt_keys = parts[8].split(":")
        # Use last sample column as tumor (convention: normal, tumor)
        tumor_idx = -1
        if len(parts) > 10:
            tumor_values = parts[10].split(":")
        else:
            tumor_values = parts[9].split(":")

        fmt_dict = dict(zip(fmt_keys, tumor_values))
        if "AF" in fmt_dict:
            try:
                vaf = float(fmt_dict["AF"].split(",")[0])
            except (ValueError, IndexError):
                pass
        if "DP" in fmt_dict:
            try:
                depth = int(fmt_dict["DP"])
            except ValueError:
                pass
        if "AD" in fmt_dict:
            try:
                ad_values = fmt_dict["AD"].split(",")
                if len(ad_values) >= 2:
                    alt_depth = int(ad_values[1])
                    if depth == 0:
                        depth = sum(int(x) for x in ad_values)
                    if vaf == 0.0 and depth > 0:
                        vaf = alt_depth / depth
            except (ValueError, IndexError):
                pass

    return SomaticVariant(
        chrom=chrom,
        pos=int(pos_str),
        ref=ref,
        alt=alt,
        gene=ann.get("gene", ""),
        transcript_id=ann.get("transcript_id", ""),
        effect=ann.get("effect", ""),
        hgvs_p=ann.get("hgvs_p", ""),
        hgvs_c=ann.get("hgvs_c", ""),
        vaf=vaf,
        depth=depth,
        alt_depth=alt_depth,
        filter_status=filt,
    )


def load_vcf(path: str | Path) -> list[SomaticVariant]:
    """Load somatic variants from an annotated VCF file.

    Args:
        path: Path to VCF file (plain text or .vcf).

    Returns:
        List of parsed SomaticVariant objects.
    """
    path = Path(path)
    variants: list[SomaticVariant] = []

    with open(path) as f:
        variants = _load_vcf_from_handle(f)

    logger.info("Loaded %d variants from %s", len(variants), path.name)
    return variants


def _load_vcf_from_handle(handle: TextIO) -> list[SomaticVariant]:
    """Parse VCF from file handle."""
    variants: list[SomaticVariant] = []
    for line in handle:
        if line.startswith("#"):
            continue
        variant = _parse_vcf_line(line)
        if variant is not None:
            variants.append(variant)
    return variants


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

# Standard coding effect types for neoantigen analysis
CODING_EFFECTS = frozenset({
    "missense_variant",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "stop_gained",
    "protein_altering_variant",
})


def filter_variants(
    variants: list[SomaticVariant],
    min_vaf: float = 0.05,
    min_depth: int = 10,
    effect_types: frozenset[str] | None = None,
    pass_only: bool = True,
) -> list[SomaticVariant]:
    """Filter somatic variants by quality and effect type.

    Args:
        variants: Input variant list.
        min_vaf: Minimum variant allele frequency.
        min_depth: Minimum read depth.
        effect_types: Allowed variant effect types. Defaults to CODING_EFFECTS.
        pass_only: Only keep variants with PASS filter status.

    Returns:
        Filtered list of variants.
    """
    if effect_types is None:
        effect_types = CODING_EFFECTS

    filtered = []
    for v in variants:
        if pass_only and v.filter_status not in ("PASS", "."):
            continue
        if v.vaf < min_vaf:
            continue
        if v.depth < min_depth:
            continue
        if effect_types and v.effect not in effect_types:
            continue
        filtered.append(v)

    logger.info(
        "Filtered %d → %d variants (min_vaf=%.2f, min_depth=%d)",
        len(variants), len(filtered), min_vaf, min_depth,
    )
    return filtered


def merge_callers(
    variant_sets: dict[str, list[SomaticVariant]],
    min_callers: int = 1,
) -> list[SomaticVariant]:
    """Merge variants from multiple callers, tracking which callers found each.

    Args:
        variant_sets: Dict mapping caller name → variant list.
        min_callers: Minimum number of callers required.

    Returns:
        Merged and deduplicated variant list.
    """
    merged: dict[str, SomaticVariant] = {}

    for caller_name, variants in variant_sets.items():
        for v in variants:
            key = v.variant_id
            if key in merged:
                merged[key].callers.add(caller_name)
                # Keep higher VAF
                if v.vaf > merged[key].vaf:
                    merged[key].vaf = v.vaf
            else:
                v.callers = {caller_name}
                merged[key] = v

    result = [v for v in merged.values() if len(v.callers) >= min_callers]
    logger.info(
        "Merged %d caller sets → %d variants (min_callers=%d)",
        len(variant_sets), len(result), min_callers,
    )
    return result
