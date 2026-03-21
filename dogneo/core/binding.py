"""MHC/DLA binding prediction wrappers.

Provides a unified interface to multiple binding prediction tools:
- NetMHCpan (MHC-I)
- NetMHCIIpan (MHC-II)
- MHCflurry (MHC-I, pan-allele)
- IEDB API (web-based, both classes)
"""

from __future__ import annotations

import csv
import json
import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import requests

from dogneo.core.peptides import MutantPeptide

logger = logging.getLogger(__name__)


@dataclass
class BindingPrediction:
    """MHC/DLA binding prediction result.

    Attributes:
        peptide_sequence: The peptide sequence evaluated.
        allele: MHC/DLA allele name.
        affinity_nm: Predicted binding affinity (nM). Lower = stronger binding.
        percentile_rank: Percentile rank (%). Lower = stronger binding.
        tool: Prediction tool used.
        mhc_class: MHC class (1 or 2).
        core_sequence: Binding core sequence (for MHC-II).
    """

    peptide_sequence: str
    allele: str
    affinity_nm: float
    percentile_rank: float = 100.0
    tool: str = "unknown"
    mhc_class: int = 1
    core_sequence: str = ""

    @property
    def is_strong_binder(self) -> bool:
        """Whether this is a strong binder (< 50 nM or < 0.5% rank)."""
        return self.affinity_nm < 50.0 or self.percentile_rank < 0.5

    @property
    def is_weak_binder(self) -> bool:
        """Whether this is at least a weak binder (< 500 nM or < 2% rank)."""
        return self.affinity_nm < 500.0 or self.percentile_rank < 2.0


# ---------------------------------------------------------------------------
# NetMHCpan Wrapper
# ---------------------------------------------------------------------------

def _write_peptide_file(peptides: list[str], path: Path) -> None:
    """Write peptide sequences to a file, one per line."""
    with open(path, "w") as f:
        for pep in peptides:
            f.write(pep + "\n")


def predict_netmhcpan(
    peptides: list[MutantPeptide],
    alleles: list[str],
    netmhcpan_path: str = "netMHCpan",
) -> list[BindingPrediction]:
    """Run NetMHCpan-4.x for MHC-I binding prediction.

    Args:
        peptides: List of mutant peptides.
        alleles: DLA/MHC allele names (e.g. ["DLA-88*001:01"]).
        netmhcpan_path: Path to netMHCpan executable.

    Returns:
        List of binding predictions.
    """
    if not peptides or not alleles:
        return []

    predictions: list[BindingPrediction] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        pep_file = Path(tmpdir) / "peptides.txt"
        _write_peptide_file([p.mut_sequence for p in peptides], pep_file)

        allele_str = ",".join(alleles)

        try:
            result = subprocess.run(
                [netmhcpan_path, "-p", str(pep_file), "-a", allele_str,
                 "-BA", "-xls", "-xlsfile", str(Path(tmpdir) / "output.xls")],
                capture_output=True, text=True, timeout=600,
            )

            output_file = Path(tmpdir) / "output.xls"
            if output_file.exists():
                predictions = _parse_netmhcpan_xls(output_file)

        except FileNotFoundError:
            logger.warning("NetMHCpan not found at: %s", netmhcpan_path)
        except subprocess.TimeoutExpired:
            logger.error("NetMHCpan timed out")
        except subprocess.CalledProcessError as e:
            logger.error("NetMHCpan error: %s", e.stderr)

    return predictions


def _parse_netmhcpan_xls(path: Path) -> list[BindingPrediction]:
    """Parse NetMHCpan XLS output format."""
    predictions: list[BindingPrediction] = []

    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header: list[str] = []
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if "Peptide" in row:
                header = row
                continue
            if not header or len(row) < len(header):
                continue

            row_dict = dict(zip(header, row))
            try:
                predictions.append(BindingPrediction(
                    peptide_sequence=row_dict.get("Peptide", ""),
                    allele=row_dict.get("MHC", row_dict.get("HLA", "")),
                    affinity_nm=float(row_dict.get("nM", row_dict.get("Aff(nM)", "99999"))),
                    percentile_rank=float(row_dict.get("%Rank", row_dict.get("Rank", "100"))),
                    tool="netmhcpan",
                    mhc_class=1,
                ))
            except (ValueError, KeyError) as e:
                logger.debug("Skipping row: %s", e)

    return predictions


# ---------------------------------------------------------------------------
# MHCflurry Wrapper
# ---------------------------------------------------------------------------

def predict_mhcflurry(
    peptides: list[MutantPeptide],
    alleles: list[str],
    mhcflurry_path: str = "mhcflurry-predict",
) -> list[BindingPrediction]:
    """Run MHCflurry for MHC-I binding prediction.

    Args:
        peptides: List of mutant peptides.
        alleles: DLA/MHC allele names.
        mhcflurry_path: Path to mhcflurry-predict executable.

    Returns:
        List of binding predictions.
    """
    if not peptides or not alleles:
        return []

    predictions: list[BindingPrediction] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        output_csv = Path(tmpdir) / "output.csv"

        # Write input CSV
        with open(input_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["allele", "peptide"])
            for allele in alleles:
                for pep in peptides:
                    writer.writerow([allele, pep.mut_sequence])

        try:
            subprocess.run(
                [mhcflurry_path, str(input_csv), "--out", str(output_csv)],
                capture_output=True, text=True, timeout=600, check=True,
            )

            with open(output_csv) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    predictions.append(BindingPrediction(
                        peptide_sequence=row["peptide"],
                        allele=row["allele"],
                        affinity_nm=float(row.get("mhcflurry_affinity", "99999")),
                        percentile_rank=float(row.get("mhcflurry_affinity_percentile", "100")),
                        tool="mhcflurry",
                        mhc_class=1,
                    ))

        except FileNotFoundError:
            logger.warning("MHCflurry not found at: %s", mhcflurry_path)
        except subprocess.CalledProcessError as e:
            logger.error("MHCflurry error: %s", e.stderr)

    return predictions


# ---------------------------------------------------------------------------
# IEDB API Wrapper
# ---------------------------------------------------------------------------

IEDB_MHC_I_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_MHC_II_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"


def predict_iedb(
    peptides: list[MutantPeptide],
    alleles: list[str],
    method: str = "recommended",
    mhc_class: int = 1,
) -> list[BindingPrediction]:
    """Query IEDB web API for MHC binding prediction.

    Note: IEDB API has rate limits. Use for small batches only.

    Args:
        peptides: List of mutant peptides.
        alleles: MHC allele names (IEDB format, e.g. "HLA-A*02:01").
        method: Prediction method (default "recommended").
        mhc_class: 1 for MHC-I, 2 for MHC-II.

    Returns:
        List of binding predictions.
    """
    if not peptides or not alleles:
        return []

    url = IEDB_MHC_I_URL if mhc_class == 1 else IEDB_MHC_II_URL
    predictions: list[BindingPrediction] = []

    sequences = "\n".join(p.mut_sequence for p in peptides)
    allele_str = ",".join(alleles)
    length_set = set(len(p.mut_sequence) for p in peptides)
    length_str = ",".join(str(l) for l in sorted(length_set))

    try:
        response = requests.post(url, data={
            "method": method,
            "sequence_text": sequences,
            "allele": allele_str,
            "length": length_str,
        }, timeout=120)

        if response.ok:
            lines = response.text.strip().split("\n")
            if len(lines) > 1:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    fields = line.split("\t")
                    if len(fields) < len(header):
                        continue
                    row = dict(zip(header, fields))
                    try:
                        predictions.append(BindingPrediction(
                            peptide_sequence=row.get("peptide", ""),
                            allele=row.get("allele", ""),
                            affinity_nm=float(row.get("ic50", "99999")),
                            percentile_rank=float(row.get("percentile_rank", "100")),
                            tool="iedb",
                            mhc_class=mhc_class,
                        ))
                    except (ValueError, KeyError):
                        continue
        else:
            logger.error("IEDB API error: %d %s", response.status_code, response.text[:200])

    except requests.RequestException as e:
        logger.error("IEDB API request failed: %s", e)

    return predictions


# ---------------------------------------------------------------------------
# Unified Interface
# ---------------------------------------------------------------------------

def predict_binding(
    peptides: list[MutantPeptide],
    alleles: list[str],
    tool: str = "netmhcpan",
    **kwargs,
) -> list[BindingPrediction]:
    """Predict MHC-I binding using specified tool.

    Args:
        peptides: Candidate mutant peptides.
        alleles: DLA/MHC allele names.
        tool: Prediction tool ("netmhcpan", "mhcflurry", "iedb").
        **kwargs: Additional tool-specific arguments.

    Returns:
        List of binding predictions.

    Raises:
        ValueError: If tool is not recognized.
    """
    if tool == "netmhcpan":
        return predict_netmhcpan(peptides, alleles, **kwargs)
    elif tool == "mhcflurry":
        return predict_mhcflurry(peptides, alleles, **kwargs)
    elif tool == "iedb":
        return predict_iedb(peptides, alleles, mhc_class=1, **kwargs)
    else:
        raise ValueError(f"Unknown binding tool: {tool}")


def predict_binding_mhcii(
    peptides: list[MutantPeptide],
    alleles: list[str],
    tool: str = "iedb",
    **kwargs,
) -> list[BindingPrediction]:
    """Predict MHC-II binding.

    Args:
        peptides: Candidate mutant peptides (15-17aa).
        alleles: DLA/MHC-II allele names.
        tool: Prediction tool ("iedb" or "netmhciipan").
        **kwargs: Additional tool-specific arguments.

    Returns:
        List of binding predictions.
    """
    if tool == "iedb":
        return predict_iedb(peptides, alleles, mhc_class=2, **kwargs)
    else:
        raise ValueError(f"MHC-II binding not supported with tool: {tool}")
