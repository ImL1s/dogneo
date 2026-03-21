"""Canine DLA (Dog Leukocyte Antigen) genotyping utilities.

Supports:
- Loading known DLA alleles from IPD-MHC database
- Running KPR tool for DLA-I genotyping from RNA-seq
- Manual allele specification
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger(__name__)


@dataclass
class DLAAllele:
    """A canine DLA (MHC) allele.

    Attributes:
        gene: DLA gene name (e.g., DLA-88, DLA-12, DLA-64).
        allele_name: Full allele name (e.g., DLA-88*001:01).
        sequence: Amino acid sequence of the allele.
        source: Data source (ipd, kpr, manual).
    """

    gene: str
    allele_name: str
    sequence: str = ""
    source: str = "manual"

    @property
    def short_name(self) -> str:
        """Short allele name without gene prefix."""
        if "*" in self.allele_name:
            return self.allele_name.split("*", 1)[1]
        return self.allele_name


# ---------------------------------------------------------------------------
# IPD-MHC Database Loading
# ---------------------------------------------------------------------------

def load_ipd_alleles(db_path: str | Path) -> list[DLAAllele]:
    """Load canine DLA alleles from IPD-MHC FASTA file.

    The IPD-MHC database provides reference sequences for all known
    DLA alleles. Download from: https://www.ebi.ac.uk/ipd/mhc/

    Expected FASTA header format:
        >DLA-88*001:01 DLA class I antigen

    Args:
        db_path: Path to IPD-MHC FASTA file for canine DLA.

    Returns:
        List of DLAAllele objects.
    """
    db_path = Path(db_path)
    alleles: list[DLAAllele] = []

    for record in SeqIO.parse(str(db_path), "fasta"):
        name = record.id
        seq = str(record.seq)

        # Parse gene name from allele name
        gene = "unknown"
        if "-" in name and "*" in name:
            gene = name.split("*")[0]

        alleles.append(DLAAllele(
            gene=gene,
            allele_name=name,
            sequence=seq,
            source="ipd",
        ))

    logger.info("Loaded %d DLA alleles from IPD-MHC: %s", len(alleles), db_path.name)
    return alleles


def filter_alleles_by_gene(
    alleles: list[DLAAllele],
    genes: list[str] | None = None,
) -> list[DLAAllele]:
    """Filter alleles by gene name.

    Args:
        alleles: Full allele list.
        genes: Gene names to keep (e.g., ["DLA-88", "DLA-12"]).

    Returns:
        Filtered allele list.
    """
    if genes is None:
        return alleles
    gene_set = set(genes)
    return [a for a in alleles if a.gene in gene_set]


# ---------------------------------------------------------------------------
# KPR Genotyping (from RNA-seq)
# ---------------------------------------------------------------------------

def run_kpr(
    rnaseq_bam: str | Path,
    reference_dir: str | Path,
    output_dir: str | Path,
    kpr_path: str = "KPR",
    threads: int = 4,
) -> list[DLAAllele]:
    """Run KPR tool to genotype DLA-I from RNA-seq data.

    KPR (K-mer based Pattern Recognition) assembles DLA-I gene sequences
    from RNA-seq reads and infers allele types.

    Reference: Zhao et al. (2023) — GitHub: VBDOL/KPR

    Args:
        rnaseq_bam: Path to RNA-seq BAM file.
        reference_dir: KPR reference database directory.
        output_dir: Output directory for KPR results.
        kpr_path: Path to KPR executable.
        threads: Number of threads.

    Returns:
        List of inferred DLAAllele objects.
    """
    rnaseq_bam = Path(rnaseq_bam)
    reference_dir = Path(reference_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = subprocess.run(
            [kpr_path,
             "--bam", str(rnaseq_bam),
             "--ref", str(reference_dir),
             "--out", str(output_dir),
             "--threads", str(threads)],
            capture_output=True, text=True, timeout=3600,
        )

        if result.returncode != 0:
            logger.error("KPR failed: %s", result.stderr)
            return []

        # Parse KPR output
        return _parse_kpr_output(output_dir)

    except FileNotFoundError:
        logger.warning("KPR not found at: %s", kpr_path)
        return []
    except subprocess.TimeoutExpired:
        logger.error("KPR timed out after 3600s")
        return []


def _parse_kpr_output(output_dir: Path) -> list[DLAAllele]:
    """Parse KPR output files for inferred DLA-I alleles.

    KPR produces a results summary with allele assignments.
    """
    alleles: list[DLAAllele] = []

    # Look for KPR result files
    result_files = list(output_dir.glob("*results*")) + list(output_dir.glob("*summary*"))

    for result_file in result_files:
        try:
            with open(result_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        allele_name = parts[0]
                        gene = allele_name.split("*")[0] if "*" in allele_name else "DLA-88"
                        alleles.append(DLAAllele(
                            gene=gene,
                            allele_name=allele_name,
                            source="kpr",
                        ))
        except OSError as e:
            logger.debug("Could not read KPR output %s: %s", result_file, e)

    logger.info("KPR inferred %d DLA-I alleles", len(alleles))
    return alleles


# ---------------------------------------------------------------------------
# Convenience Functions
# ---------------------------------------------------------------------------

def infer_dla_from_rnaseq(
    fastq_paths: list[str | Path],
    reference_dir: str | Path = "",
    output_dir: str | Path = "/tmp/kpr_output",
) -> list[str]:
    """Convenience function to infer DLA alleles from RNA-seq FASTQ files.

    Args:
        fastq_paths: RNA-seq FASTQ files (will be aligned first).
        reference_dir: KPR reference directory.
        output_dir: Output directory.

    Returns:
        List of allele name strings.
    """
    logger.warning(
        "infer_dla_from_rnaseq: Full implementation requires BAM alignment first. "
        "Please provide pre-aligned BAM to run_kpr() directly."
    )
    return []


def parse_allele_string(allele_str: str) -> list[DLAAllele]:
    """Parse comma-separated allele string into DLAAllele objects.

    Args:
        allele_str: Comma-separated allele names (e.g., "DLA-88*001:01,DLA-88*501:01").

    Returns:
        List of DLAAllele objects.
    """
    alleles: list[DLAAllele] = []
    for name in allele_str.split(","):
        name = name.strip()
        if not name:
            continue
        gene = name.split("*")[0] if "*" in name else "DLA-88"
        alleles.append(DLAAllele(
            gene=gene,
            allele_name=name,
            source="manual",
        ))
    return alleles
