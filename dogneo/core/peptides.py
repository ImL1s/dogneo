"""Mutant peptide window generation from somatic variants.

Generates candidate neoantigen peptide sequences centered on mutation sites,
for both MHC class I (8-11aa) and MHC class II (15-17aa) binding prediction.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from dogneo.core.variants import SomaticVariant

logger = logging.getLogger(__name__)


@dataclass
class MutantPeptide:
    """A candidate neoantigen peptide with wild-type and mutant sequences.

    Attributes:
        gene: Gene symbol.
        variant_id: Source variant identifier.
        mutation: Protein-level mutation description.
        wt_sequence: Wild-type peptide sequence.
        mut_sequence: Mutant peptide sequence.
        position: Mutation position within the peptide (0-indexed).
        length: Peptide length.
        mhc_class: MHC class (1 or 2).
    """

    gene: str
    variant_id: str
    mutation: str
    wt_sequence: str
    mut_sequence: str
    position: int
    length: int
    mhc_class: int = 1

    @property
    def peptide_id(self) -> str:
        """Unique identifier for this peptide."""
        return f"{self.variant_id}|{self.mut_sequence}|{self.length}mer"

    def is_self_different(self) -> bool:
        """Whether mutant differs from wild-type."""
        return self.wt_sequence != self.mut_sequence


class ProteinDatabase:
    """Simple protein sequence database loaded from FASTA.

    Maps transcript IDs or gene names to protein sequences for
    peptide window extraction.
    """

    def __init__(self) -> None:
        self._by_transcript: dict[str, str] = {}
        self._by_gene: dict[str, str] = {}

    def load_fasta(self, path: str | Path) -> None:
        """Load protein sequences from FASTA file.

        Supported header formats:
        - Ensembl: >TRANSCRIPT gene_symbol:SYMBOL gene:ENSG... (uses gene_symbol)
        - Simple:  >TRANSCRIPT gene=SYMBOL or gene:SYMBOL

        When multiple isoforms exist per gene, keeps the longest sequence
        to maximize peptide coverage.

        Args:
            path: Path to protein FASTA file.
        """
        path = Path(path)
        count = 0
        for record in SeqIO.parse(str(path), "fasta"):
            transcript_id = record.id
            seq_str = str(record.seq)
            self._by_transcript[transcript_id] = seq_str

            # Try to extract gene name from description
            # Supports: gene_symbol:SYMBOL (Ensembl, preferred),
            #           gene=SYMBOL, gene:SYMBOL (simple formats)
            desc = record.description
            gene_name = None
            gene_fallback = None
            for part in desc.split():
                if part.startswith("gene_symbol:"):
                    gene_name = part.split(":", 1)[1]
                    break  # gene_symbol is authoritative, stop
                if part.startswith("gene="):
                    gene_fallback = part.split("=", 1)[1]
                elif part.startswith("gene:") and not part.startswith("gene_biotype:"):
                    val = part.split(":", 1)[1]
                    # Skip Ensembl gene IDs (e.g. ENSCAFG00000016714.4)
                    if not val.startswith("ENS"):
                        gene_fallback = val
            gene_name = gene_name or gene_fallback
            if gene_name:
                # Keep longest isoform per gene for best peptide coverage
                if gene_name not in self._by_gene or len(seq_str) > len(self._by_gene[gene_name]):
                    self._by_gene[gene_name] = seq_str
            count += 1

        logger.info("Loaded %d protein sequences from %s", count, path.name)

    def get_sequence(self, transcript_id: str = "", gene: str = "") -> str | None:
        """Look up protein sequence by transcript ID or gene name.

        Args:
            transcript_id: Transcript identifier.
            gene: Gene symbol (fallback).

        Returns:
            Protein sequence string, or None if not found.
        """
        if transcript_id and transcript_id in self._by_transcript:
            return self._by_transcript[transcript_id]
        if gene and gene in self._by_gene:
            return self._by_gene[gene]
        return None


def _parse_missense_hgvsp(hgvs_p: str) -> tuple[str, int, str] | None:
    """Parse a missense HGVS protein notation.

    Examples:
        p.V600E -> ('V', 600, 'E')
        p.Ala123Thr -> ('A', 123, 'T')

    Returns:
        Tuple of (ref_aa, position, alt_aa) or None if unparseable.
    """
    # Standard three-letter to one-letter mapping
    aa3to1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*",
    }

    p = hgvs_p.replace("p.", "").strip()
    if not p:
        return None

    # Try single-letter format: V600E
    if len(p) >= 3 and p[0].isalpha() and p[-1].isalpha():
        ref_aa = p[0]
        alt_aa = p[-1]
        pos_str = p[1:-1]
        try:
            pos = int(pos_str)
            return ref_aa, pos, alt_aa
        except ValueError:
            pass

    # Try three-letter format: Ala123Thr
    for ref_3, ref_1 in aa3to1.items():
        if p.startswith(ref_3):
            rest = p[len(ref_3):]
            for alt_3, alt_1 in aa3to1.items():
                if rest.endswith(alt_3):
                    pos_str = rest[:-len(alt_3)]
                    try:
                        pos = int(pos_str)
                        return ref_1, pos, alt_1
                    except ValueError:
                        pass
            break

    return None


def _generate_windows(
    protein_seq: str,
    mut_pos: int,
    ref_aa: str,
    alt_aa: str,
    lengths: list[int],
) -> list[tuple[str, str, int]]:
    """Generate all peptide windows of given lengths centered on mutation.

    Args:
        protein_seq: Full protein sequence.
        mut_pos: 1-based mutation position in protein.
        ref_aa: Reference amino acid.
        alt_aa: Alternate amino acid.
        lengths: Peptide lengths to generate.

    Returns:
        List of (wt_peptide, mut_peptide, mut_position_in_peptide).
    """
    idx = mut_pos - 1  # Convert to 0-based
    if idx < 0 or idx >= len(protein_seq):
        return []

    results: list[tuple[str, str, int]] = []

    for length in lengths:
        # Slide window: mutation can be at any position in the peptide
        for offset in range(length):
            start = idx - offset
            end = start + length

            if start < 0 or end > len(protein_seq):
                continue

            wt_peptide = protein_seq[start:end]
            # Create mutant peptide
            mut_peptide_list = list(wt_peptide)
            mut_peptide_list[offset] = alt_aa
            mut_peptide = "".join(mut_peptide_list)

            # Skip if identical (synonymous at peptide level)
            if wt_peptide == mut_peptide:
                continue

            results.append((wt_peptide, mut_peptide, offset))

    return results


def generate_peptides(
    variant: SomaticVariant,
    protein_db: ProteinDatabase,
    lengths: list[int] | None = None,
) -> list[MutantPeptide]:
    """Generate candidate MHC-I neoantigen peptides for a missense variant.

    Args:
        variant: Annotated somatic variant with HGVS protein notation.
        protein_db: Protein sequence database.
        lengths: Peptide lengths (default [8, 9, 10, 11] for MHC-I).

    Returns:
        List of MutantPeptide candidates.
    """
    if lengths is None:
        lengths = [8, 9, 10, 11]

    if not variant.hgvs_p:
        return []

    parsed = _parse_missense_hgvsp(variant.hgvs_p)
    if parsed is None:
        logger.debug("Could not parse HGVS.p: %s", variant.hgvs_p)
        return []

    ref_aa, pos, alt_aa = parsed

    protein_seq = protein_db.get_sequence(
        transcript_id=variant.transcript_id,
        gene=variant.gene,
    )
    if protein_seq is None:
        logger.debug("No protein sequence found for %s / %s",
                      variant.transcript_id, variant.gene)
        return []

    windows = _generate_windows(protein_seq, pos, ref_aa, alt_aa, lengths)

    peptides = []
    for wt_seq, mut_seq, mut_pos in windows:
        peptides.append(MutantPeptide(
            gene=variant.gene,
            variant_id=variant.variant_id,
            mutation=variant.hgvs_p,
            wt_sequence=wt_seq,
            mut_sequence=mut_seq,
            position=mut_pos,
            length=len(mut_seq),
            mhc_class=1,
        ))

    return peptides


def generate_mhcii_peptides(
    variant: SomaticVariant,
    protein_db: ProteinDatabase,
    lengths: list[int] | None = None,
) -> list[MutantPeptide]:
    """Generate candidate MHC-II neoantigen peptides (longer windows).

    Args:
        variant: Annotated somatic variant.
        protein_db: Protein sequence database.
        lengths: Peptide lengths (default [15, 16, 17] for MHC-II).

    Returns:
        List of MutantPeptide candidates.
    """
    if lengths is None:
        lengths = [15, 16, 17]

    peptides = generate_peptides(variant, protein_db, lengths)
    for p in peptides:
        p.mhc_class = 2
    return peptides
