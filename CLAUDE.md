# Dog mRNA SOS — Canine Neoantigen Prioritization Tool

## Overview
Research-use-only (RUO) pipeline for canine tumor neoantigen prioritization.
Chains: tumor/normal sequencing → somatic variants → peptide generation → DLA binding → immunogenicity ranking.

## Tech Stack
- **Language**: Python 3.10+
- **Workflow**: Snakemake
- **Reference Genome**: CanFam3.1 (default)
- **License**: Apache-2.0

## Project Structure
```
dogneo/
├── core/          # VCF parsing, peptides, binding, ranking, DLA typing
├── pipeline/      # Snakemake workflow definitions
├── llm/           # Pluggable LLM layer (CLI / local / cloud)
├── report/        # HTML/Markdown report generation
├── export/        # FASTA/TSV/JSON output
└── cli.py         # Click-based CLI entrypoint
```

## Commands
```bash
# Install
pip install -e ".[all]"

# Run tests
pytest tests/ -v

# Lint
ruff check dogneo/ tests/

# CLI
dogneo --help
dogneo run --config config.yaml
dogneo rank --vcf somatic.vcf --alleles DLA-88*001:01
dogneo report --input ranked.json --output report.html
```

## Conventions
- All outputs must carry RUO disclaimers
- Type hints required on all public functions
- Docstrings in Google style
- Tests in `tests/` mirroring `dogneo/` structure
- LLM backends: prefer CLI (free) → local → cloud API

## LLM Backend Priority
1. **CLI** (Gemini CLI, Claude Code, Codex) — free, no API key
2. **Local** (llama-cpp-python GGUF) — offline, privacy
3. **Cloud API** (OpenAI, Anthropic, Gemini) — highest quality
