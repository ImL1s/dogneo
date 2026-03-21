# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

DogNeo is a research-use-only (RUO) pipeline for canine tumor neoantigen prioritization. It chains: annotated VCF → somatic variant parsing → mutant peptide generation → DLA binding prediction → multi-factor immunogenicity ranking → mRNA construct design → multi-format export + interactive UI.

**All outputs must carry the RUO disclaimer** (defined in `dogneo/__init__.py:RUO_DISCLAIMER`).

## Commands

```bash
# Install (dev)
pip install -e ".[all]"

# Run all tests
pytest tests/ -v

# Run a single test file / single test
pytest tests/test_variants.py -v
pytest tests/test_variants.py::test_load_vcf -v

# Lint
ruff check dogneo/ tests/

# Type check
mypy dogneo/

# CLI (10 commands)
dogneo setup                  # download CanFam3.1 proteome (~15 MB, one-time)
dogneo demo                   # full pipeline on bundled osteosarcoma data
dogneo rank --vcf somatic.vcf # rank with auto IEDB binding prediction
dogneo rerank --candidates results/candidates.json --binding netmhcpan.tsv
dogneo design-mrna --candidates results/candidates.json --top-n 10
dogneo report --input ranked.json --output report.html
dogneo ui                     # launch Streamlit interactive dashboard
dogneo ui --port 9000 --demo  # custom port, pre-load demo
dogneo check-llm              # show available LLM backends
dogneo run --config pipeline.yaml  # full Snakemake pipeline
```

## Architecture

```
dogneo/
├── app/                    # Application service layer
│   ├── rank_pipeline.py    #   RankInput → run_rank_pipeline() → RankResult
│   └── rerank_pipeline.py  #   Import external binding + re-score
├── core/                   # Computational modules
│   ├── variants.py         #   VCF parsing (SnpEff ANN / VEP CSQ)
│   ├── peptides.py         #   Mutant peptide window generation
│   ├── binding.py          #   NetMHCpan / MHCflurry wrappers
│   ├── iedb_client.py      #   IEDB API client (cache + rate limit)
│   ├── ranking.py          #   6-factor composite scoring
│   ├── mrna_designer.py    #   Canine codon optimization + mRNA construct
│   ├── expression.py       #   RNA-seq expression (Salmon/Kallisto)
│   └── dla_typing.py       #   DLA genotyping
├── ui/                     # Streamlit frontend
│   ├── app.py              #   4-page dashboard (Upload/Ranking/Viz/Report)
│   └── charts.py           #   Plotly chart builders (bar/heatmap/radar/scatter)
├── llm/                    # Pluggable LLM layer
│   ├── router.py           #   3-tier routing (CLI → Local → Cloud)
│   ├── explainer.py        #   Step-by-step plain-language explanations
│   ├── backends.py         #   Backend implementations
│   ├── cli_wrapper.py      #   Subprocess wrapper for AI CLIs
│   └── prompts.py          #   Prompt templates
├── pipeline/               # Snakemake workflow
├── report/                 # HTML/Markdown generation
├── export/                 # FASTA/TSV/JSON exporters
├── data/                   # Reference + demo data
└── cli.py                  # Click CLI (thin — calls app/ services)
```

### Key Architecture Decisions

- **`app/` service layer** — `run_rank_pipeline()` is the single entry point for ranking, called by CLI, UI, and tests. CLI is thin (parse args → call service → print output).
- **Binding fallback chain** — `auto` resolves: netMHCpan (local) → IEDB API (remote, free) → none (offline). Implemented in `_resolve_binding_tool()`.
- **IEDB client** — `core/iedb_client.py` caches results to `~/.dogneo/cache/iedb/` with rate limiting (1s between requests). Supports DLA alleles via NetMHCpan EL backend.
- **Explainer hooks** — When `--llm-tier` is set, `PipelineExplainer` generates plain-language explanations at each pipeline step, stored in `RankResult.explanations`.

## Data Flow

```
SomaticVariant (variants.py)
  → parsed from VCF (SnpEff ANN or VEP CSQ)
  → filtered by VAF ≥ 0.05, depth ≥ 10, coding effects

MutantPeptide (peptides.py)
  → sliding window over CanFam3.1 proteome
  → MHC-I: 8-11aa; MHC-II: 15-17aa

BindingPrediction (iedb_client.py or binding.py)
  → IEDB API (default, free) or NetMHCpan/MHCflurry (local)
  → <50 nM = strong binder, <500 nM = weak binder

NeoantigenCandidate (ranking.py)
  → 6-factor composite score → ranked list

MRNAConstruct (mrna_designer.py)
  → top candidates → codon-optimized cassette + UTR + poly-A
  → FASTA output for synthesis
```

## Scoring System

`ranking.py:rank_candidates()` computes a weighted composite score (0–1):

| Factor | Weight | Logic |
|--------|--------|-------|
| binding_affinity | 0.30 | Log-linear: 50nM→1.0, 500nM→0.5, 5000nM→0.0 |
| expression | 0.20 | Log-linear on TPM |
| vaf | 0.15 | Linear: VAF×2, capped at 1.0 |
| self_difference | 0.15 | AA differences between WT/mut peptide |
| agretopicity | 0.10 | WT/mut binding ratio (tanh-normalized) |
| caller_agreement | 0.10 | num_callers / 3 |

## LLM Layer

Three-tier fallback in `llm/router.py`. LLMs are used for:
1. **Report summaries** — narrative interpretation in HTML reports
2. **Step explanations** — `llm/explainer.py` explains each pipeline step in plain language
3. **mRNA design review** — commentary on construct design (core optimization is algorithmic)

All computational analysis works without any LLM. `--llm-tier none` (default) skips LLM entirely.

## Key Domain Concepts

- **DLA** (Dog Leukocyte Antigen): canine MHC. 6 DLA-88 Class I alleles bundled from IPD-MHC.
- **HGVS.p**: protein mutation notation (e.g., `p.Val245Ile`). Parsed in `peptides.py`.
- **CanFam3.1**: default canine proteome, auto-downloaded via `dogneo setup`.
- **IEDB API**: free web service for DLA binding prediction. Used as default when NetMHCpan unavailable.
- **Canine codon table**: `mrna_designer.py:CANINE_CODON_TABLE` — optimal codons for Canis lupus familiaris.

## Test Structure

189 tests in `tests/` covering: variants, peptides, ranking, exporters, LLM router, CLI, rank pipeline, IEDB client, binding integration, rerank pipeline, explainer, mRNA designer, UI charts.

Shared fixtures in `conftest.py`. Test data in `tests/data/`; demo data in `dogneo/data/demo/`.

## Conventions

- Python 3.10+, type hints on all public functions, Google-style docstrings
- Ruff: line-length 100, target py310, rules E/F/W/I/N/UP/B/SIM
- CI runs on Python 3.11 and 3.12
- `build/` directory is a stale setuptools copy — source of truth is `dogneo/`
- Optional deps: `[llm]`, `[bio]`, `[pipeline]`, `[ui]`, `[dev]`, `[all]`
