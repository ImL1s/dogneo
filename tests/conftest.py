"""Shared pytest fixtures for DogNeo test suite."""
from __future__ import annotations

from pathlib import Path

import pytest

from dogneo.core.binding import BindingPrediction
from dogneo.core.peptides import MutantPeptide
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.core.variants import SomaticVariant


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def sample_vcf_path() -> Path:
    """Path to synthetic canine VCF with SnpEff/VEP annotations."""
    return DATA_DIR / "sample.vcf"


@pytest.fixture
def sample_expression_path() -> Path:
    """Path to synthetic Salmon quant.sf file."""
    return DATA_DIR / "salmon_quant.sf"


@pytest.fixture
def sample_alleles() -> list[str]:
    """DLA-88 allele names from IPD-MHC database."""
    path = DATA_DIR / "alleles.txt"
    return [line.strip() for line in path.read_text().strip().splitlines()]


@pytest.fixture
def sample_variant() -> SomaticVariant:
    """A single TP53 missense variant for unit tests."""
    return SomaticVariant(
        chrom="chr1",
        pos=15000100,
        ref="G",
        alt="A",
        gene="TP53",
        transcript_id="ENSCAFT00000013077",
        effect="missense_variant",
        hgvs_p="p.Val245Ile",
        hgvs_c="c.733G>A",
        vaf=0.40,
        depth=50,
        alt_depth=20,
        filter_status="PASS",
        expression_tpm=85.3,
    )


@pytest.fixture
def sample_variants(sample_variant: SomaticVariant) -> list[SomaticVariant]:
    """Three test variants: TP53 (missense), BRAF (missense), synonymous CDH1."""
    braf = SomaticVariant(
        chrom="chr5",
        pos=32500200,
        ref="T",
        alt="C",
        gene="BRAF",
        transcript_id="ENSCAFT00000005678",
        effect="missense_variant",
        hgvs_p="p.Val600Ala",
        hgvs_c="c.1799T>C",
        vaf=0.375,
        depth=40,
        alt_depth=15,
        filter_status="PASS",
        expression_tpm=2.1,
    )
    cdh1 = SomaticVariant(
        chrom="chr2",
        pos=22000400,
        ref="A",
        alt="G",
        gene="CDH1",
        transcript_id="ENSCAFT00000006666",
        effect="synonymous_variant",
        hgvs_p="p.Arg200=",
        hgvs_c="c.600A>G",
        vaf=0.375,
        depth=40,
        alt_depth=15,
        filter_status="PASS",
    )
    return [sample_variant, braf, cdh1]


@pytest.fixture
def sample_peptide() -> MutantPeptide:
    """9-mer mutant peptide from TP53 p.Val245Ile."""
    return MutantPeptide(
        gene="TP53",
        variant_id="chr1:15000100:G>A",
        mutation="p.Val245Ile",
        wt_sequence="RAVVGAPPS",
        mut_sequence="RAIVGAPPS",
        position=2,
        length=9,
        mhc_class=1,
    )


@pytest.fixture
def sample_binding() -> BindingPrediction:
    """Strong binding prediction for TP53 peptide."""
    return BindingPrediction(
        peptide_sequence="RAIVGAPPS",
        allele="DLA-88*001:01",
        affinity_nm=42.5,
        percentile_rank=0.8,
        tool="NetMHCpan",
        mhc_class=1,
    )


@pytest.fixture
def sample_candidate(
    sample_variant: SomaticVariant,
    sample_peptide: MutantPeptide,
    sample_binding: BindingPrediction,
) -> NeoantigenCandidate:
    """Complete NeoantigenCandidate for TP53 V245I."""
    return NeoantigenCandidate(
        variant=sample_variant,
        peptide=sample_peptide,
        binding=sample_binding,
        expression_tpm=85.3,
        composite_score=0.85,
        rank=1,
        score_components={
            "binding_affinity": 0.95,
            "expression": 0.88,
            "vaf": 0.80,
            "self_difference": 0.33,
            "agretopicity": 0.50,
            "caller_agreement": 0.33,
        },
    )


@pytest.fixture
def tmp_output_dir(tmp_path: Path) -> Path:
    """Temporary output directory, pre-created."""
    out = tmp_path / "output"
    out.mkdir()
    return out
