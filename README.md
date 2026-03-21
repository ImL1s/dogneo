<div align="center">

# 🧬 DogNeo

**Canine Neoantigen Prioritization Pipeline**

From somatic variants to ranked vaccine candidates — purpose-built for dogs.

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/tests-137_passed-brightgreen.svg)](#-verification)
[![install](https://img.shields.io/badge/install-GitHub-orange.svg)](#-installation)

[Quick Start](#-quick-start) · [Usage](#-usage) · [LLM Config](#-llm-backend-configuration) · [Architecture](#-pipeline-architecture) · [License](#-license)

</div>

---

> **⚠️ FOR RESEARCH USE ONLY (RUO)** — This tool is intended solely for computational research purposes.
> It does not provide clinical diagnoses, veterinary advice, or treatment recommendations.

## Why DogNeo?

Cancer immunotherapy has shown remarkable results in human medicine — but dogs get cancer too. When a team [used AI + mRNA to treat a dog's tumor](https://people.com/owner-uses-ai-to-create-vaccine-for-dogs-cancer-11683328), it highlighted a gap: there was **no open-source, canine-specific neoantigen pipeline**.

DogNeo fills that gap. It takes annotated somatic variants from a dog's tumor, generates candidate mutant peptides, predicts binding to **DLA (Dog Leukocyte Antigen)** molecules, and ranks neoantigens — giving researchers a starting point for personalized canine cancer vaccines.

```
Tumor/Normal DNA ──→ Variant Calling ──→ Peptide Generation ──→ DLA Binding ──→ Ranked Candidates
       RNA-seq ──→ Expression Quantification ─────────────────────┘
       DLA typing (KPR) ──────────────────────────────────────────┘
```

## ✨ Key Features

| Feature | Description |
|---------|-------------|
| 🐕 **Canine-specific** | DLA alleles from [IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/), CanFam3.1 proteome, KPR genotyping |
| 🔬 **Full pipeline** | BWA → GATK → Mutect2/Strelka → VEP/SnpEff → Salmon → NetMHCpan → Ranking |
| 🧪 **Out-of-the-box demo** | 3 commands to run a complete analysis on real canine osteosarcoma data |
| 🤖 **AI-assisted reports** | 3-tier LLM backend (CLI free → Local → Cloud) for narrative interpretation |
| 📦 **Multi-format output** | TSV, JSON, FASTA, and interactive HTML reports |
| ♻️ **Reproducible** | Snakemake workflow with Docker containerization |

## 🚀 Quick Start

```bash
pip install git+https://github.com/ImL1s/dogneo.git   # 1. Install
dogneo setup                                           # 2. Download reference data (~15 MB)
dogneo demo                                            # 3. Run demo pipeline ✨
```

That's it. The demo runs a full pipeline on **bundled canine osteosarcoma data** (8 published mutations across TP53, BRAF, KRAS, PIK3CA, PTEN) and generates ranked candidates in TSV, JSON, and FASTA formats.

<details>
<summary><strong>📋 Expected output</strong></summary>

```
🧬 DogNeo v0.1.0 — Running demo pipeline
📂 Using bundled demo data
📂 Using reference: ~/.dogneo/data/CanFam3.1.pep.all.fa
   8 total → 6 coding variants
   228 peptides from 6 variants
   912 unscored candidates

✅ Demo complete! Results: dogneo_demo_results/CANINE_OSA_DEMO
   📄 candidates.tsv   — Tab-separated candidates
   📄 candidates.json  — Structured JSON with metadata
   📄 candidates.fasta — Peptide sequences for wet lab
```

</details>

## 📦 Installation

### Prerequisites

| Tool | Version | Check | Notes |
|------|---------|-------|-------|
| **Python** | ≥ 3.10 | `python --version` | Required |
| **pip** | latest | `pip --version` | Included with Python |
| **git** | any | `git --version` | Required |

### Install from GitHub (Recommended)

```bash
pip install git+https://github.com/ImL1s/dogneo.git
```

### With extras

```bash
pip install "dogneo[bio] @ git+https://github.com/ImL1s/dogneo.git"     # + pysam, pyvcf3
pip install "dogneo[llm] @ git+https://github.com/ImL1s/dogneo.git"     # + openai, anthropic
pip install "dogneo[all] @ git+https://github.com/ImL1s/dogneo.git"     # Everything
```

### From source (development)

```bash
git clone https://github.com/ImL1s/dogneo.git
cd dogneo
pip install -e ".[all]"
```

### Reference Data (~15 MB, one-time)

```bash
dogneo setup
```

This downloads the CanFam3.1 proteome from Ensembl FTP to `~/.dogneo/data/`. DLA alleles are bundled with the package — no additional download needed.

## 📖 Usage

### Rank neoantigens from your own VCF

```bash
# Minimal — auto-detects cached proteome + bundled DLA alleles
dogneo rank --vcf my_somatic_variants.vcf

# Full options
dogneo rank \
  --vcf somatic_annotated.vcf \
  --expression salmon_quant.sf \
  --alleles "DLA-88*001:01,DLA-88*501:01" \
  --protein-db CanFam3.1.pep.all.fa \
  --output-dir results/ \
  --formats tsv,json,fasta,html \
  --llm-tier cli
```

### Generate AI-assisted report

```bash
# From existing candidates JSON
dogneo report \
  --input results/candidates.json \
  --output report.html \
  --llm-tier cli    # free, uses Gemini CLI / Claude Code
```

### Run full Snakemake pipeline

```bash
dogneo run --config pipeline_config.yaml
```

### CLI Reference

| Command | Description |
|---------|-------------|
| `dogneo setup` | Download reference data (CanFam3.1 proteome) |
| `dogneo demo` | Run full pipeline on bundled demo data |
| `dogneo rank` | Rank neoantigens from a VCF file |
| `dogneo report` | Generate HTML/Markdown report |
| `dogneo check-llm` | Display status of all LLM backends |
| `dogneo version` | Show version |

## 🤖 LLM Backend Configuration

The AI report layer is **entirely optional** — all computational analysis works without any LLM. LLMs are only used for generating narrative interpretations in HTML reports.

### Three tiers, prioritized by cost

| Tier | Backend | Cost | Privacy | Setup |
|------|---------|------|---------|-------|
| 1. CLI | Gemini CLI, Claude Code, Codex | **Free** | Data sent to cloud | Install any AI CLI tool |
| 2. Local | llama-cpp (GGUF models) | **Free** | **Fully offline** | Download a GGUF model |
| 3. Cloud | OpenAI, Anthropic, Google Gemini | Pay per token | Data sent to cloud | Set API key |

### Check available backends

```bash
$ dogneo check-llm

📡 CLI Backends (Tier 1 — Free):
   ✅ gemini
   ✅ claude
   ✅ codex
💾 Local Backends (Tier 2 — Offline):
   ⚪ No local model configured
☁️  Cloud Backends (Tier 3 — API):
   ❌ OpenAI  (OPENAI_API_KEY)
   ❌ Anthropic  (ANTHROPIC_API_KEY)
```

### Environment Variables (Cloud tier only)

```bash
# Optional — only needed if using cloud LLM tier
export OPENAI_API_KEY=sk-...
export ANTHROPIC_API_KEY=sk-ant-...
export GOOGLE_API_KEY=AI...
```

## 🔬 Pipeline Architecture

```
                          ┌──────────────┐
                          │  Tumor DNA   │
                          │  Normal DNA  │
                          │  RNA-seq     │
                          └──────┬───────┘
                                 │
                    ┌────────────┼────────────┐
                    ▼            ▼             ▼
             ┌───────────┐ ┌──────────┐ ┌──────────┐
             │ BWA-MEM   │ │   STAR   │ │  Salmon  │
             │ alignment │ │ alignment│ │ quant.   │
             └─────┬─────┘ └────┬─────┘ └────┬─────┘
                   │            │             │
             ┌─────┴─────┐     │        Expression
             │   GATK    │     │          TPMs
             │ Mutect2   │     │             │
             │ Strelka2  │     │             │
             └─────┬─────┘     │             │
                   │           │             │
             ┌─────┴─────┐    │             │
             │ VEP/SnpEff│    │             │
             │ annotation│    │             │
             └─────┬─────┘    │             │
                   │          │             │
                   ▼          ▼             ▼
            ┌──────────────────────────────────┐
            │       DogNeo Core Engine         │
            │  ┌─────────────────────────────┐ │
            │  │ Variant Parsing & Filtering │ │
            │  │ Mutant Peptide Generation   │ │
            │  │ DLA Binding Prediction      │ │
            │  │ Multi-factor Ranking        │ │
            │  └─────────────────────────────┘ │
            └──────────────┬───────────────────┘
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
         ┌────────┐  ┌────────┐  ┌──────────┐
         │  TSV   │  │  JSON  │  │  HTML    │
         │  FASTA │  │        │  │ + LLM   │
         └────────┘  └────────┘  └──────────┘
```

### Computational Steps

| Step | Tool | Description |
|------|------|-------------|
| 1 | BWA-MEM / STAR | Sequence alignment (DNA / RNA) |
| 2 | GATK | MarkDuplicates, BQSR |
| 3 | Mutect2 + Strelka2 | Somatic variant calling |
| 4 | VEP / SnpEff | Variant annotation |
| 5 | Salmon / Kallisto | Expression quantification |
| 6 | DogNeo | Mutant peptide windows (8–11aa MHC-I, 15–17aa MHC-II) |
| 7 | KPR | DLA genotyping from RNA-seq |
| 8 | NetMHCpan / MHCflurry / IEDB | Binding prediction |
| 9 | DogNeo | Multi-factor immunogenicity scoring |
| 10 | DogNeo | TSV / JSON / FASTA / HTML report |

## 🗂️ Project Structure

```
dogneo/
├── core/                 # Core computational modules
│   ├── variants.py       #   VCF parsing & somatic variant handling
│   ├── peptides.py       #   Mutant peptide window generation
│   ├── expression.py     #   RNA-seq expression quantification
│   ├── binding.py        #   MHC/DLA binding prediction wrappers
│   ├── ranking.py        #   Multi-factor neoantigen scoring
│   └── dla_typing.py     #   Canine DLA genotyping
├── data/                 # Reference & demo data
│   ├── manager.py        #   Auto-download CanFam3.1 from Ensembl FTP
│   └── demo/             #   Bundled demo data (VCF, expression, alleles)
├── pipeline/             # Snakemake workflow definitions
│   ├── Snakefile         #   Main workflow
│   └── config.yaml       #   Pipeline configuration template
├── llm/                  # Pluggable LLM layer (optional)
│   ├── router.py         #   3-tier routing (CLI → Local → Cloud)
│   ├── backends.py       #   OpenAI / Anthropic / Gemini backends
│   ├── cli_wrapper.py    #   Subprocess wrapper for AI CLIs
│   └── prompts.py        #   Prompt templates for neoantigen analysis
├── report/               # Report generation
│   └── generator.py      #   HTML/Markdown report with AI summary
├── export/               # Output format converters
│   └── exporters.py      #   FASTA, TSV, JSON exporters
└── cli.py                # Click-based CLI entrypoint
```

## 📊 Reference Data

| Resource | Source | Notes |
|----------|--------|-------|
| Canine proteome | [Ensembl](https://ftp.ensembl.org/) CanFam3.1 | Auto-downloaded via `dogneo setup` (~15 MB) |
| DLA alleles | [IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/) | 6 DLA-88 Class I alleles, bundled |
| DLA typing | KPR ([paper](https://doi.org/10.1016/j.isci.2023.105996)) | RNA-seq based MHC-I genotyping |
| Canine genome | CanFam3.1, GSD_1.0, UMICH_Zoey_3.1 | Multiple assemblies supported |

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

```bash
# Development setup
git clone https://github.com/ImL1s/dogneo.git
cd dogneo
pip install -e ".[all]"
python -m pytest tests/ -v        # Run tests (137 passing)
```

## ⚖️ Legal Notices

<details>
<summary><strong>RUO Disclaimer</strong></summary>

All outputs are labeled **"FOR RESEARCH USE ONLY — NOT FOR CLINICAL OR VETERINARY DIAGNOSTIC USE."**
This tool does not provide medical advice, diagnoses, or treatment recommendations.
Users must obtain appropriate ethical approvals before using this tool with real animal samples.

</details>

<details>
<summary><strong>Patent Notice</strong></summary>

Existing patents (e.g., ASU WO2018223094A1) cover canine cancer vaccine methods.
This tool performs **computational analysis only** and does not constitute vaccine design or manufacture.

</details>

## 🙏 Acknowledgments

DogNeo stands on the shoulders of incredible open-source tools:

**Upstream:** [BWA](https://github.com/lh3/bwa) · [GATK](https://gatk.broadinstitute.org/) · [Mutect2](https://gatk.broadinstitute.org/) · [Strelka2](https://github.com/Illumina/strelka) · [STAR](https://github.com/alexdobin/STAR) · [Salmon](https://combine-lab.github.io/salmon/) · [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) · [MHCflurry](https://github.com/openvax/mhcflurry) · [VEP](https://www.ensembl.org/vep) · [SnpEff](https://pcingola.github.io/SnpEff/)

**Inspired by:** [pVACtools](https://github.com/griffithlab/pVACtools) · [OpenVax/Vaxrank](https://github.com/openvax/vaxrank) · [nextNEOpi](https://github.com/icbi-lab/nextNEOpi) · [MiroFish](https://github.com/666ghj/MiroFish) · [Ollama](https://github.com/ollama/ollama)

**Motivated by:** The story of [Paul Conyngham using AI to develop an mRNA treatment](https://people.com/owner-uses-ai-to-create-vaccine-for-dogs-cancer-11683328) for his dog Rosie's cancer — and the realization that accessible tools can accelerate veterinary oncology research.

## 📄 License

[Apache License 2.0](LICENSE)

---

<p align="center">
  <sub>Made with 🐕 for dogs everywhere fighting cancer.</sub>
</p>
