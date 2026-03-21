"""Global configuration for DogNeo pipeline."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class ReferenceConfig:
    """Reference genome and annotation configuration."""

    genome_fasta: str = ""
    genome_name: str = "CanFam3.1"
    annotation_gtf: str = ""
    known_sites_vcf: str = ""
    star_index_dir: str = ""
    bwa_index_prefix: str = ""
    snpeff_db: str = "CanFam3.1.105"


@dataclass
class ToolConfig:
    """External bioinformatics tool paths."""

    bwa: str = "bwa"
    star: str = "STAR"
    gatk: str = "gatk"
    samtools: str = "samtools"
    salmon: str = "salmon"
    snpeff: str = "snpEff"
    vep: str = "vep"
    netmhcpan: str = "netMHCpan"
    netmhciipan: str = "netMHCIIpan"
    mhcflurry: str = "mhcflurry-predict"


@dataclass
class LLMConfig:
    """LLM backend configuration."""

    # Tier priority: cli -> local -> cloud
    default_tier: str = "cli"

    # CLI backends
    gemini_cli_model: str = "gemini-2.5-flash"
    claude_cli_model: str = "claude-sonnet-4-6"
    codex_cli_model: str = "gpt-5.4"
    cli_timeout: int = 180

    # Local backend
    local_model_path: str = ""
    local_n_ctx: int = 4096
    local_n_gpu_layers: int = -1

    # Cloud API keys (read from environment)
    openai_api_key: str = ""
    anthropic_api_key: str = ""
    google_api_key: str = ""
    cloud_model: str = "gpt-4o"

    def __post_init__(self) -> None:
        if not self.openai_api_key:
            self.openai_api_key = os.environ.get("OPENAI_API_KEY", "")
        if not self.anthropic_api_key:
            self.anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
        if not self.google_api_key:
            self.google_api_key = os.environ.get("GOOGLE_API_KEY", "")


@dataclass
class PipelineConfig:
    """Main pipeline configuration."""

    sample_id: str = ""
    tumor_dna_fastq: list[str] = field(default_factory=list)
    normal_dna_fastq: list[str] = field(default_factory=list)
    tumor_rna_fastq: list[str] = field(default_factory=list)
    dla_alleles: list[str] = field(default_factory=list)
    output_dir: str = "results"
    threads: int = 4

    # Peptide generation
    mhci_peptide_lengths: list[int] = field(default_factory=lambda: [8, 9, 10, 11])
    mhcii_peptide_lengths: list[int] = field(default_factory=lambda: [15, 16, 17])

    # Filtering
    min_vaf: float = 0.05
    min_expression_tpm: float = 1.0
    max_binding_affinity_nm: float = 500.0

    # Sub-configs
    reference: ReferenceConfig = field(default_factory=ReferenceConfig)
    tools: ToolConfig = field(default_factory=ToolConfig)
    llm: LLMConfig = field(default_factory=LLMConfig)


def load_config(path: str | Path) -> PipelineConfig:
    """Load pipeline configuration from YAML file.

    Args:
        path: Path to YAML configuration file.

    Returns:
        Populated PipelineConfig instance.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path) as f:
        raw: dict[str, Any] = yaml.safe_load(f) or {}

    config = PipelineConfig()

    # Top-level fields
    for key in ("sample_id", "output_dir", "threads", "min_vaf",
                "min_expression_tpm", "max_binding_affinity_nm"):
        if key in raw:
            setattr(config, key, raw[key])

    # List fields
    for key in ("tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq",
                "dla_alleles", "mhci_peptide_lengths", "mhcii_peptide_lengths"):
        if key in raw:
            setattr(config, key, raw[key])

    # Sub-configs
    if "reference" in raw:
        for k, v in raw["reference"].items():
            if hasattr(config.reference, k):
                setattr(config.reference, k, v)

    if "tools" in raw:
        for k, v in raw["tools"].items():
            if hasattr(config.tools, k):
                setattr(config.tools, k, v)

    if "llm" in raw:
        for k, v in raw["llm"].items():
            if hasattr(config.llm, k):
                setattr(config.llm, k, v)

    return config
