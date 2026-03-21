# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

DogNeo is a research-use-only (RUO) pipeline for canine tumor neoantigen prioritization. It chains: annotated VCF → somatic variant parsing → mutant peptide generation → DLA binding prediction → multi-factor immunogenicity ranking → multi-format export (TSV/JSON/FASTA/HTML).

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

# CLI
dogneo --help
dogneo setup                  # download CanFam3.1 proteome (~15 MB, one-time)
dogneo demo                   # full pipeline on bundled osteosarcoma data
dogneo rank --vcf somatic.vcf --alleles "DLA-88*001:01"
dogneo report --input ranked.json --output report.html
dogneo check-llm              # show available LLM backends
```

## Architecture — Data Flow

The core pipeline is a linear chain of dataclasses, each building on the previous:

```
SomaticVariant (variants.py)
  → parsed from VCF (SnpEff ANN or VEP CSQ fields)
  → filtered by VAF, depth, coding effect type

MutantPeptide (peptides.py)
  → sliding window over ProteinDatabase sequences
  → requires HGVS.p (missense only) + protein FASTA from CanFam3.1
  → MHC-I: 8-11aa windows; MHC-II: 15-17aa windows

BindingPrediction (binding.py)
  → wraps NetMHCpan / MHCflurry / IEDB API
  → key thresholds: <50 nM = strong binder, <500 nM = weak binder

NeoantigenCandidate (ranking.py)
  → combines variant + peptide + binding + expression
  → composite_score from 6 weighted factors (see Scoring below)
  → exported via exporters.py (FASTA/TSV/JSON) and generator.py (HTML/Markdown)
```

The CLI (`cli.py`) orchestrates this chain. The `rank` command is the main entry point for standalone analysis; `demo` invokes `rank` programmatically via Click's `CliRunner`.

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

Weights are in `ScoringWeights` dataclass — adjustable per invocation.

## LLM Layer

Three-tier fallback system in `llm/router.py`:
1. **CLI** — subprocess calls to gemini/claude/codex CLIs (free, no API key)
2. **Local** — llama-cpp-python with GGUF models (offline)
3. **Cloud** — OpenAI/Anthropic/Gemini APIs (requires env vars)

LLMs are **only** used for narrative HTML report summaries. All computational analysis works without any LLM. The `--llm-tier none` default skips LLM entirely.

## Key Domain Concepts

- **DLA** (Dog Leukocyte Antigen): canine equivalent of human HLA. Alleles bundled in `dogneo/data/demo/alleles.txt` (6 DLA-88 Class I alleles from IPD-MHC).
- **HGVS.p notation**: protein-level mutation format (e.g., `p.Val245Ile`). Parsed by `peptides.py:_parse_missense_hgvsp()` supporting both single-letter (V245I) and three-letter (Val245Ile) formats.
- **CanFam3.1**: default canine reference proteome, auto-downloaded to `~/.dogneo/data/` via `dogneo setup`.
- **VCF annotation**: supports both SnpEff (ANN field) and VEP (CSQ field) formats. Parsing in `variants.py`.

## Test Structure

Tests in `tests/` mirror `dogneo/core/` modules. Shared fixtures in `conftest.py` provide:
- `sample_variant` / `sample_variants` — synthetic canine variants (TP53, BRAF, CDH1)
- `sample_peptide`, `sample_binding`, `sample_candidate` — complete test objects
- `sample_vcf_path`, `sample_expression_path` — paths under `tests/data/`

Test data files live in `tests/data/`; demo data (bundled with package) lives in `dogneo/data/demo/`.

## Conventions

- Python 3.10+, type hints on all public functions, Google-style docstrings
- Ruff config: line-length 100, target py310, rules: E/F/W/I/N/UP/B/SIM
- CI runs on Python 3.11 and 3.12 (see `.github/workflows/ci.yml`)
- `build/` directory contains a stale setuptools copy — source of truth is `dogneo/`

## LLM Backend Priority

1. **CLI** (Gemini CLI, Claude Code, Codex) — free, no API key
2. **Local** (llama-cpp-python GGUF) — offline, privacy
3. **Cloud API** (OpenAI, Anthropic, Gemini) — highest quality
