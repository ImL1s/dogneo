"""Tests for dogneo.core.variants — VCF loading and variant filtering."""
from __future__ import annotations

from pathlib import Path

import pytest

from dogneo.core.variants import (
    SomaticVariant,
    load_vcf,
    filter_variants,
    merge_callers,
    _parse_info_field,
    _extract_snpeff_annotation,
    _extract_vep_annotation,
    CODING_EFFECTS,
)


# ---------------------------------------------------------------------------
# INFO / ANN / CSQ parsing
# ---------------------------------------------------------------------------

class TestInfoParsing:
    """Unit tests for VCF INFO field parsing helpers."""

    def test_parse_info_kv(self):
        result = _parse_info_field("DP=50;MQ=60;ANN=allele_ann")
        assert result["DP"] == "50"
        assert result["MQ"] == "60"
        assert result["ANN"] == "allele_ann"

    def test_parse_info_flag(self):
        result = _parse_info_field("SOMATIC;DP=30")
        assert result["SOMATIC"] == "true"
        assert result["DP"] == "30"


class TestSnpEffAnnotation:
    """Tests for SnpEff ANN field extraction."""

    def test_full_ann(self):
        ann = "A|missense_variant|MODERATE|TP53|ENSG1|transcript|ENST1|coding|7/11|c.733G>A|p.Val245Ile|1001|733|245"
        result = _extract_snpeff_annotation(ann)
        assert result["effect"] == "missense_variant"
        assert result["gene"] == "TP53"
        assert result["transcript_id"] == "ENST1"
        assert result["hgvs_c"] == "c.733G>A"
        assert result["hgvs_p"] == "p.Val245Ile"

    def test_short_ann_returns_empty(self):
        ann = "A|missense|MODERATE|TP53"  # only 4 fields
        result = _extract_snpeff_annotation(ann)
        assert result == {}


class TestVepAnnotation:
    """Tests for VEP CSQ field extraction."""

    def test_full_csq(self):
        csq = "T|missense_variant|MODERATE|KRAS||transcript|ENST7||||c.182A>T|p.Gln61Leu"
        result = _extract_vep_annotation(csq)
        assert result["effect"] == "missense_variant"
        assert result["gene"] == "KRAS"
        assert result["hgvs_c"] == "c.182A>T"
        assert result["hgvs_p"] == "p.Gln61Leu"

    def test_short_csq_returns_empty(self):
        csq = "T|missense|MODERATE"
        result = _extract_vep_annotation(csq)
        assert result == {}

    def test_csq_missing_hgvs_fields(self):
        """VEP with only 8 fields → hgvs_c / hgvs_p should be empty."""
        csq = "T|missense_variant|MODERATE|KRAS||transcript|ENST7|coding"
        result = _extract_vep_annotation(csq)
        assert result["gene"] == "KRAS"
        assert result["hgvs_c"] == ""
        assert result["hgvs_p"] == ""


# ---------------------------------------------------------------------------
# VCF loading from file
# ---------------------------------------------------------------------------

class TestLoadVcf:
    """Tests for loading variants from synthetic VCF."""

    def test_load_snpeff_variants(self, sample_vcf_path: Path):
        variants = load_vcf(sample_vcf_path)
        # VCF has 8 data lines total
        assert len(variants) == 8

    def test_snpeff_tp53_parsed(self, sample_vcf_path: Path):
        variants = load_vcf(sample_vcf_path)
        tp53 = [v for v in variants if v.gene == "TP53"]
        assert len(tp53) == 1
        assert tp53[0].effect == "missense_variant"
        assert tp53[0].hgvs_p == "p.Val245Ile"
        assert tp53[0].hgvs_c == "c.733G>A"

    def test_vep_kras_parsed(self, sample_vcf_path: Path):
        """VEP CSQ entry should also be parsed (KRAS line)."""
        variants = load_vcf(sample_vcf_path)
        kras = [v for v in variants if v.gene == "KRAS"]
        assert len(kras) == 1
        assert kras[0].hgvs_p == "p.Gln61Leu"

    def test_vaf_extraction(self, sample_vcf_path: Path):
        variants = load_vcf(sample_vcf_path)
        tp53 = [v for v in variants if v.gene == "TP53"][0]
        assert tp53.vaf == pytest.approx(0.40, abs=0.01)
        assert tp53.depth == 50
        assert tp53.alt_depth == 20

    def test_frameshift_detected(self, sample_vcf_path: Path):
        variants = load_vcf(sample_vcf_path)
        pten = [v for v in variants if v.gene == "PTEN"]
        assert len(pten) == 1
        assert pten[0].effect == "frameshift_variant"

    def test_synonymous_detected(self, sample_vcf_path: Path):
        variants = load_vcf(sample_vcf_path)
        cdh1 = [v for v in variants if v.gene == "CDH1"]
        assert len(cdh1) == 1
        assert cdh1[0].effect == "synonymous_variant"


# ---------------------------------------------------------------------------
# Variant ID format
# ---------------------------------------------------------------------------

class TestVariantId:
    """SomaticVariant.variant_id formatting."""

    def test_format(self, sample_variant: SomaticVariant):
        assert sample_variant.variant_id == "chr1:15000100:G>A"

    def test_is_coding_missense(self, sample_variant: SomaticVariant):
        assert sample_variant.is_coding is True

    def test_is_coding_synonymous(self):
        v = SomaticVariant(chrom="chr2", pos=100, ref="A", alt="G", effect="synonymous_variant")
        assert v.is_coding is False


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

class TestFilterVariants:
    """Tests for variant filtering logic."""

    def test_filter_by_vaf(self, sample_variants: list[SomaticVariant]):
        # All have VAF >= 0.375 → all pass with min_vaf=0.3
        result = filter_variants(sample_variants, min_vaf=0.3, min_depth=0, pass_only=False)
        # Synonymous is excluded by effect filter (CODING_EFFECTS)
        assert all(v.effect in CODING_EFFECTS for v in result)

    def test_filter_excludes_low_vaf(self):
        low_vaf = SomaticVariant(
            chrom="chr1", pos=100, ref="A", alt="T",
            effect="missense_variant", vaf=0.02, depth=50,
            filter_status="PASS",
        )
        result = filter_variants([low_vaf], min_vaf=0.05, min_depth=0)
        assert len(result) == 0

    def test_filter_excludes_synonymous(self, sample_variants: list[SomaticVariant]):
        result = filter_variants(sample_variants, min_vaf=0.0, min_depth=0, pass_only=False)
        genes = [v.gene for v in result]
        assert "CDH1" not in genes

    def test_pass_only(self):
        """LowQual-filtered variants should be excluded by default."""
        v = SomaticVariant(
            chrom="chr4", pos=500, ref="G", alt="C",
            effect="missense_variant", vaf=0.25, depth=20,
            filter_status="LowQual",
        )
        result = filter_variants([v], min_vaf=0.0, min_depth=0, pass_only=True)
        assert len(result) == 0

    def test_pass_dot_allowed(self):
        """'.' filter status should be treated as PASS."""
        v = SomaticVariant(
            chrom="chr4", pos=500, ref="G", alt="C",
            effect="missense_variant", vaf=0.25, depth=20,
            filter_status=".",
        )
        result = filter_variants([v], min_vaf=0.0, min_depth=0, pass_only=True)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Merge callers
# ---------------------------------------------------------------------------

class TestMergeCallers:
    """Tests for multi-caller variant merging."""

    def test_single_caller(self, sample_variant: SomaticVariant):
        result = merge_callers({"mutect2": [sample_variant]})
        assert len(result) == 1
        assert "mutect2" in result[0].callers

    def test_two_callers(self, sample_variant: SomaticVariant):
        v2 = SomaticVariant(
            chrom=sample_variant.chrom, pos=sample_variant.pos,
            ref=sample_variant.ref, alt=sample_variant.alt,
            gene="TP53", vaf=0.38, depth=55,
        )
        result = merge_callers({"mutect2": [sample_variant], "strelka": [v2]})
        assert len(result) == 1
        assert result[0].callers == {"mutect2", "strelka"}
        # Should keep higher VAF
        assert result[0].vaf == pytest.approx(0.40, abs=0.01)

    def test_min_callers_filter(self, sample_variant: SomaticVariant):
        result = merge_callers({"mutect2": [sample_variant]}, min_callers=2)
        assert len(result) == 0
