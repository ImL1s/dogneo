# 🧬 Dog mRNA SOS — Canine Neoantigen Prioritization Tool

> **⚠️ RESEARCH USE ONLY (RUO)** — This tool is intended solely for computational
> research purposes. It does not provide clinical diagnoses, veterinary advice, or
> treatment recommendations. Not a medical device. Not for clinical use.

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)

An open-source pipeline for canine tumor **neoantigen prioritization**, linking
tumor/normal sequencing data to candidate neoantigen generation and ranking.

## Overview

This tool automates the computational workflow for identifying potential
neoantigens from canine tumor samples:

```
Tumor/Normal DNA → Variant Calling → Peptide Generation → DLA Binding → Ranked Candidates
         RNA-seq → Expression Quantification ─────────────────┘
         DLA typing (KPR) ──────────────────────────────────────┘
```

## Key Features

- **Full pipeline**: BWA → GATK → Mutect2/Strelka → VEP/SnpEff → Salmon → NetMHCpan → Ranking
- **Canine-specific**: DLA (Dog Leukocyte Antigen) typing via KPR, IPD-MHC database support
- **Pluggable LLM**: 3-tier AI backend (CLI free → Local → Cloud API) for report generation
- **Reproducible**: Snakemake workflow with Docker containerization
- **Multi-format output**: TSV, JSON, FASTA with provenance tracking

## Installation

```bash
# Basic install (includes BioPython, Click, Pandas, requests)
pip install -e .

# With bioinformatics extras (pysam, pyvcf3)
pip install -e ".[bio]"

# With LLM support
pip install -e ".[llm]"

# Everything
pip install -e ".[all]"
```

## Quick Start (3 commands)

```bash
# 1. Install
pip install -e .

# 2. Download reference data (one-time, ~15 MB)
dogneo setup

# 3. Run demo — zero config needed
dogneo demo
```

That's it! The demo runs a full pipeline on bundled canine osteosarcoma data
(8 published mutations from TP53, BRAF, KRAS, PIK3CA, PTEN) and generates
TSV, JSON, and FASTA output.

### Advanced Usage

```bash
# Rank from your own VCF (auto-detects cached proteome + bundled DLA alleles)
dogneo rank --vcf my_somatic_variants.vcf

# With specific alleles and binding prediction
dogneo rank \
  --vcf somatic_annotated.vcf \
  --expression salmon_quant.sf \
  --alleles "DLA-88*001:01,DLA-88*501:01" \
  --output-dir results

# Run full Snakemake pipeline
dogneo run --config pipeline_config.yaml

# Generate HTML report from candidates
dogneo report --input candidates.json --output report.html
```

## LLM Backend Configuration

Three tiers of AI backends, prioritized by cost:

| Tier | Backend | API Key? | Use Case |
|------|---------|----------|----------|
| 1. CLI | Gemini CLI, Claude Code, Codex | ❌ Free | Default for all tasks |
| 2. Local | llama-cpp (GGUF) | ❌ Offline | Privacy-sensitive data |
| 3. Cloud | OpenAI, Anthropic, Gemini API | ✅ Required | Highest quality reports |

```bash
# Use CLI backend (default, free)
dogneo report --model gemini-cli --input candidates.json

# Use local model
dogneo report --model local:mistral-7b.gguf --input candidates.json

# Use cloud API
dogneo report --model openai:gpt-4o --input candidates.json
```

## Project Structure

```
dogneo/
├── core/               # Core computational modules
│   ├── variants.py     # VCF parsing & somatic variant handling
│   ├── peptides.py     # Mutant peptide window generation
│   ├── expression.py   # RNA-seq expression quantification
│   ├── binding.py      # MHC/DLA binding prediction wrappers
│   ├── ranking.py      # Multi-factor neoantigen scoring
│   └── dla_typing.py   # Canine DLA genotyping
├── pipeline/           # Snakemake workflow definitions
│   ├── Snakefile       # Main workflow (10 steps)
│   └── config.yaml     # Pipeline configuration template
├── llm/                # Pluggable LLM layer
│   ├── router.py       # 3-tier routing (CLI → Local → Cloud)
│   ├── backends.py     # Backend implementations
│   ├── cli_wrapper.py  # Subprocess wrapper for AI CLIs
│   └── prompts.py      # Prompt templates
├── report/             # Report generation
│   ├── generator.py    # Jinja2-based HTML/Markdown reports
│   └── templates/      # Report templates
├── export/             # Output format converters
│   └── exporters.py    # FASTA, TSV, JSON exporters
└── cli.py              # Click-based CLI entrypoint
```

## Reference Data

- **Canine genome**: CanFam3.1 (default), GSD_1.0, UMICH_Zoey_3.1
- **DLA alleles**: [IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/) canine database
- **DLA typing**: [KPR](https://github.com/VBDOL/KPR) from RNA-seq

## Computational Steps

1. **Sequence alignment** — BWA-MEM (DNA), STAR (RNA)
2. **GATK preprocessing** — MarkDuplicates, BQSR
3. **Variant calling** — Mutect2 + Strelka2
4. **Variant annotation** — VEP / SnpEff
5. **Expression quantification** — Salmon / Kallisto
6. **Peptide generation** — Mutant peptide windows (8-11aa MHC-I, 15-17aa MHC-II)
7. **DLA typing** — KPR genotyping from RNA-seq
8. **Binding prediction** — NetMHCpan / MHCflurry / IEDB
9. **Immunogenicity ranking** — Multi-factor scoring
10. **Output** — TSV/JSON/FASTA + HTML report

## Legal Notices

> **Patent Notice**: Existing patents (e.g., ASU WO2018223094A1) cover canine
> cancer vaccine methods. This tool performs computational analysis only and does
> not constitute vaccine design or manufacture.

> **RUO Disclaimer**: All outputs are labeled "FOR RESEARCH USE ONLY — NOT FOR
> CLINICAL OR VETERINARY DIAGNOSTIC USE." Users must obtain appropriate ethical
> approvals before using this tool with real animal samples.

## Acknowledgments

Built upon open-source bioinformatics tools: BWA, GATK, Mutect2, Strelka2,
STAR, Salmon, NetMHCpan, MHCflurry, SnpEff/VEP. Inspired by human neoantigen
pipelines: pVACtools, OpenVax/Vaxrank, and the VetClaw veterinary AI skill library.

## License

Apache License 2.0 — See [LICENSE](LICENSE) for details.
