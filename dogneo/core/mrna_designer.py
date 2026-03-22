"""mRNA construct designer for canine neoantigen vaccines.

Generates codon-optimized multi-epitope mRNA constructs with UTRs,
signal peptides, linker sequences, and poly-A tails.
Uses canine (Canis lupus familiaris) codon usage table.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dogneo.core.ranking import NeoantigenCandidate

if TYPE_CHECKING:
    from dogneo.llm.router import LLMRouter

logger = logging.getLogger(__name__)

# Canine optimal codons (Canis lupus familiaris)
# Source: Codon Usage Database (Kazusa) — most frequent codon per amino acid
CANINE_CODON_TABLE: dict[str, str] = {
    "A": "GCC",  # Ala
    "R": "AGG",  # Arg
    "N": "AAC",  # Asn
    "D": "GAC",  # Asp
    "C": "UGC",  # Cys
    "E": "GAG",  # Glu
    "Q": "CAG",  # Gln
    "G": "GGC",  # Gly
    "H": "CAC",  # His
    "I": "AUC",  # Ile
    "L": "CUG",  # Leu
    "K": "AAG",  # Lys
    "M": "AUG",  # Met
    "F": "UUC",  # Phe
    "P": "CCC",  # Pro
    "S": "AGC",  # Ser
    "T": "ACC",  # Thr
    "W": "UGG",  # Trp
    "Y": "UAC",  # Tyr
    "V": "GUG",  # Val
}

# Standard RNA codon table for translation verification
_RNA_CODON_TO_AA: dict[str, str] = {
    "GCC": "A", "GCU": "A", "GCA": "A", "GCG": "A",
    "AGG": "R", "AGA": "R", "CGG": "R", "CGC": "R", "CGA": "R", "CGU": "R",
    "AAC": "N", "AAU": "N",
    "GAC": "D", "GAU": "D",
    "UGC": "C", "UGU": "C",
    "GAG": "E", "GAA": "E",
    "CAG": "Q", "CAA": "Q",
    "GGC": "G", "GGU": "G", "GGA": "G", "GGG": "G",
    "CAC": "H", "CAU": "H",
    "AUC": "I", "AUU": "I", "AUA": "I",
    "CUG": "L", "CUC": "L", "CUA": "L", "CUU": "L", "UUA": "L", "UUG": "L",
    "AAG": "K", "AAA": "K",
    "AUG": "M",
    "UUC": "F", "UUU": "F",
    "CCC": "P", "CCU": "P", "CCA": "P", "CCG": "P",
    "AGC": "S", "AGU": "S", "UCG": "S", "UCC": "S", "UCA": "S", "UCU": "S",
    "ACC": "T", "ACU": "T", "ACA": "T", "ACG": "T",
    "UGG": "W",
    "UAC": "Y", "UAU": "Y",
    "GUG": "V", "GUC": "V", "GUA": "V", "GUU": "V",
    "UAA": "*", "UAG": "*", "UGA": "*",
}

# Default UTR sequences (based on human beta-globin, widely used in mRNA vaccines)
DEFAULT_5UTR = "AGAAUAAACUAGUAUUCUUCUGGUCCCCACAGACUCAGAGAGAACCCGCCACC"
DEFAULT_3UTR = (
    "GCUGGAGCCUCGGUGGCCAUGCUUCUUGCCCCUUGGGCCUCCCCCCAGCCCCUCCUCCCCUUCCUGCACCCGUACCCCC"
    "GUGGUCUUUGAAUAAAGUCUGAGUGGGCGGC"
)

# Common linker sequences for multi-epitope vaccines
LINKER_GPGPG = "GPGPG"  # Flexible linker, preserves epitope processing
LINKER_AAY = "AAY"       # Short linker for proteasomal cleavage


def codon_optimize(protein_seq: str, species: str = "canine") -> str:
    """Reverse-translate protein to codon-optimized mRNA.

    Args:
        protein_seq: Amino acid sequence (single-letter codes).
        species: Target species for codon optimization.

    Returns:
        Codon-optimized RNA sequence (ACGU alphabet).
    """
    codons = []
    for aa in protein_seq.upper():
        if aa in CANINE_CODON_TABLE:
            codons.append(CANINE_CODON_TABLE[aa])
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    return "".join(codons)


def _translate_rna(rna_seq: str) -> str:
    """Translate RNA sequence back to protein (for verification).

    Args:
        rna_seq: RNA sequence (ACGU alphabet).

    Returns:
        Protein sequence (single-letter codes).
    """
    protein = []
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i + 3]
        aa = _RNA_CODON_TO_AA.get(codon, "?")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


@dataclass
class MRNAConstruct:
    """Complete mRNA construct ready for synthesis.

    Structure: 5'cap — 5'UTR — signal_peptide — [epitope-linker]×N — 3'UTR — poly(A)
    """

    five_prime_utr: str
    three_prime_utr: str
    signal_peptide: str = ""
    antigen_cassette: str = ""
    linker_sequences: list[str] = field(default_factory=list)
    epitope_sequences: list[str] = field(default_factory=list)
    epitope_genes: list[str] = field(default_factory=list)
    poly_a_length: int = 120
    epitope_count: int = 0

    @property
    def full_sequence(self) -> str:
        """Assemble the complete mRNA sequence."""
        parts = [self.five_prime_utr]
        if self.signal_peptide:
            parts.append(self.signal_peptide)
        parts.append(self.antigen_cassette)
        parts.append(self.three_prime_utr)
        parts.append("A" * self.poly_a_length)
        return "".join(parts)

    def to_fasta(self, name: str = "DogNeo_mRNA_construct") -> str:
        """Export as FASTA format."""
        seq = self.full_sequence
        header = (
            f">{name} "
            f"epitopes={self.epitope_count} "
            f"length={len(seq)}nt "
            f"genes={','.join(self.epitope_genes)}"
        )
        # Wrap at 80 characters
        lines = [header]
        for i in range(0, len(seq), 80):
            lines.append(seq[i:i + 80])
        return "\n".join(lines)


def design_construct(
    candidates: list[NeoantigenCandidate],
    top_n: int = 10,
    linker: str = LINKER_GPGPG,
    llm_router: LLMRouter | None = None,
) -> MRNAConstruct:
    """Design multi-epitope mRNA construct from top candidates.

    Selects top N candidates, concatenates their mutant peptide sequences
    with linker regions, and wraps in UTRs and poly-A tail.

    Args:
        candidates: Ranked neoantigen candidates.
        top_n: Number of top candidates to include.
        linker: Linker peptide sequence between epitopes.
        llm_router: Optional LLM for design review.

    Returns:
        MRNAConstruct ready for export.
    """
    selected = candidates[:top_n]

    # Deduplicate by peptide sequence
    seen_peptides: set[str] = set()
    unique_epitopes: list[tuple[str, str]] = []  # (peptide_seq, gene)
    for c in selected:
        if c.peptide.mut_sequence not in seen_peptides:
            seen_peptides.add(c.peptide.mut_sequence)
            unique_epitopes.append((c.peptide.mut_sequence, c.variant.gene))

    # Build antigen cassette: epitope1-linker-epitope2-linker-...
    cassette_protein_parts: list[str] = []
    linker_seqs: list[str] = []
    for i, (pep_seq, _gene) in enumerate(unique_epitopes):
        cassette_protein_parts.append(pep_seq)
        if i < len(unique_epitopes) - 1:
            cassette_protein_parts.append(linker)
            linker_seqs.append(linker)

    cassette_protein = "".join(cassette_protein_parts)

    # Codon-optimize the cassette
    cassette_rna = codon_optimize(cassette_protein)

    return MRNAConstruct(
        five_prime_utr=DEFAULT_5UTR,
        three_prime_utr=DEFAULT_3UTR,
        antigen_cassette=cassette_rna,
        linker_sequences=linker_seqs,
        epitope_sequences=[pep for pep, _ in unique_epitopes],
        epitope_genes=[gene for _, gene in unique_epitopes],
        poly_a_length=120,
        epitope_count=len(unique_epitopes),
    )
