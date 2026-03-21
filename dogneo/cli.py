"""DogNeo CLI — Canine Neoantigen Prioritization Tool.

FOR RESEARCH USE ONLY — NOT FOR CLINICAL OR VETERINARY DIAGNOSTIC USE.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from dogneo import __version__, RUO_DISCLAIMER

logger = logging.getLogger("dogneo")


@click.group()
@click.version_option(__version__, prog_name="dogneo")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose logging.")
def cli(verbose: bool) -> None:
    """🧬 DogNeo — Canine Neoantigen Prioritization Tool.

    FOR RESEARCH USE ONLY.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)-7s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True),
              help="Pipeline config YAML file.")
@click.option("--cores", "-j", default=4, help="Number of cores for Snakemake.")
@click.option("--dry-run", "-n", is_flag=True, help="Snakemake dry run (don't execute).")
def run(config: str, cores: int, dry_run: bool) -> None:
    """Run the full Snakemake pipeline.

    Executes all steps from alignment through neoantigen ranking.
    """
    click.echo(f"⚠️  {RUO_DISCLAIMER}")
    click.echo(f"🧬 DogNeo v{__version__} — Starting pipeline")

    snakefile = Path(__file__).parent / "pipeline" / "Snakefile"
    if not snakefile.exists():
        click.secho("❌ Snakefile not found!", fg="red")
        sys.exit(1)

    import subprocess
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", config,
        "--cores", str(cores),
    ]
    if dry_run:
        cmd.append("-n")

    click.echo(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


@cli.command()
@click.option("--vcf", required=True, type=click.Path(exists=True),
              help="Annotated VCF file (SnpEff/VEP).")
@click.option("--expression", type=click.Path(exists=True),
              help="Salmon/Kallisto quantification file.")
@click.option("--sample-id", default="SAMPLE", help="Sample identifier.")
@click.option("--output-dir", "-o", default="results", help="Output directory.")
@click.option("--alleles", default="",
              help="Comma-separated DLA alleles (e.g., DLA-88*001:01,DLA-88*002:01).")
@click.option("--mhci-lengths", default="8,9,10,11",
              help="MHC-I peptide lengths.")
@click.option("--mhcii-lengths", default="15,16,17",
              help="MHC-II peptide lengths.")
@click.option("--formats", default="tsv,json",
              help="Output formats: tsv,json,fasta,html.")
@click.option("--llm-tier", type=click.Choice(["cli", "local", "cloud", "none"]),
              default="none", help="LLM tier for AI-assisted analysis.")
def rank(
    vcf: str,
    expression: str | None,
    sample_id: str,
    output_dir: str,
    alleles: str,
    mhci_lengths: str,
    mhcii_lengths: str,
    formats: str,
    llm_tier: str,
) -> None:
    """Rank neoantigens from an annotated VCF.

    Performs peptide generation, binding prediction, scoring,
    and exports results in the specified formats.
    """
    click.echo(f"⚠️  {RUO_DISCLAIMER}")
    click.echo(f"🧬 DogNeo v{__version__} — Ranking mode")

    from dogneo.core.variants import load_vcf, filter_variants
    from dogneo.core.peptides import generate_peptides
    from dogneo.core.expression import load_expression, annotate_expression
    from dogneo.core.ranking import NeoantigenCandidate, rank_candidates
    from dogneo.export.exporters import export_tsv, export_json, export_fasta

    outdir = Path(output_dir) / sample_id
    outdir.mkdir(parents=True, exist_ok=True)
    format_list = [f.strip().lower() for f in formats.split(",")]

    # Parse alleles
    allele_list = [a.strip() for a in alleles.split(",") if a.strip()] if alleles else []

    # Step 1: Load variants
    click.echo("📂 Loading variants from VCF...")
    variants = load_vcf(vcf)
    coding = filter_variants(variants)
    click.echo(f"   {len(variants)} total → {len(coding)} coding variants")

    # Step 2: Load expression (optional)
    expr_data = None
    if expression:
        click.echo("📊 Loading expression data...")
        expr_data = load_expression(expression)
        coding = annotate_expression(coding, expr_data)

    # Step 3: Generate peptides
    click.echo("🔬 Generating mutant peptides...")
    mhci_lens = [int(x) for x in mhci_lengths.split(",")]
    click.echo(f"   Peptide lengths: MHC-I {mhci_lens}")

    # NOTE: full peptide generation requires a canine protein FASTA database.
    # In the full pipeline (Snakemake), this is handled by the alignment steps.
    from dogneo.core.peptides import ProteinDatabase, generate_peptides
    from dogneo.core.binding import BindingPrediction
    from dogneo.core.ranking import build_candidates

    protein_db = ProteinDatabase()
    # TODO: accept --protein-db CLI flag for standalone usage
    peptides_by_variant: dict[str, list] = {}
    predictions_by_peptide: dict[str, list] = {}

    for v in coding:
        peps = generate_peptides(v, protein_db, lengths=mhci_lens)
        if peps:
            peptides_by_variant[v.variant_id] = peps

    if not peptides_by_variant:
        click.secho(
            "⚠️  No peptides generated — this likely means no canine protein DB "
            "was loaded. Use the full Snakemake pipeline (dogneo run) or provide "
            "a pre-built candidates JSON to the report command.",
            fg="yellow",
        )

    # Step 4: Ranking
    click.echo("📊 Scoring and ranking candidates...")
    candidates = build_candidates(coding, peptides_by_variant, predictions_by_peptide)

    if allele_list and candidates:
        ranked = rank_candidates(candidates)
        click.echo(f"   {len(ranked)} candidates ranked")
    else:
        ranked = candidates
        click.echo(f"   {len(ranked)} candidates (unranked — no alleles or binding data)")

    # Step 5: Export
    if "tsv" in format_list:
        export_tsv(candidates, outdir / "candidates.tsv")
    if "json" in format_list:
        export_json(candidates, outdir / "candidates.json", sample_id)
    if "fasta" in format_list:
        export_fasta(candidates, outdir / "candidates.fasta")
    if "html" in format_list:
        from dogneo.report.generator import ReportGenerator
        llm_router = None
        if llm_tier != "none":
            from dogneo.config import LLMConfig
            from dogneo.llm.router import LLMRouter
            llm_config = LLMConfig(default_tier=llm_tier)
            llm_router = LLMRouter(config=llm_config)

        gen = ReportGenerator(llm_router=llm_router)
        gen.generate_html(
            candidates, sample_id,
            parameters={"VCF": vcf, "Alleles": alleles, "Formats": formats},
            alleles=allele_list,
            output_path=outdir / "report.html",
        )

    click.secho(f"✅ Results written to: {outdir}", fg="green")


@cli.command()
@click.option("--input", "-i", "input_path", required=True,
              type=click.Path(exists=True), help="Candidates JSON file.")
@click.option("--format", "-f", "fmt", default="html",
              type=click.Choice(["html", "markdown"]), help="Report format.")
@click.option("--output", "-o", required=True, help="Output file path.")
@click.option("--llm-tier", type=click.Choice(["cli", "local", "cloud", "none"]),
              default="none", help="LLM tier for AI summary.")
def report(input_path: str, fmt: str, output: str, llm_tier: str) -> None:
    """Generate a report from candidates JSON."""
    import json as _json
    click.echo(f"⚠️  {RUO_DISCLAIMER}")

    with open(input_path) as f:
        data = _json.load(f)

    total = data.get("metadata", {}).get("total_candidates", "?")
    sample_id = data.get("metadata", {}).get("sample_id", "UNKNOWN")
    click.echo(f"📄 Generating {fmt} report from {total} candidates...")

    from dogneo.report.generator import ReportGenerator
    from dogneo.config import LLMConfig
    from dogneo.llm.router import LLMRouter

    llm_router = None
    if llm_tier != "none":
        llm_config = LLMConfig(default_tier=llm_tier)
        llm_router = LLMRouter(config=llm_config)

    gen = ReportGenerator(llm_router=llm_router)
    output_path = Path(output)

    # The JSON's "candidates" list already has serialized candidate dicts
    candidate_dicts = data.get("candidates", [])

    if fmt == "html":
        gen.generate_html(
            [], sample_id,
            parameters=data.get("metadata", {}).get("parameters", {}),
            alleles=data.get("metadata", {}).get("alleles", []),
            output_path=output_path,
            pre_rendered_candidates=candidate_dicts,
        )
    else:
        gen.generate_markdown(
            [], sample_id,
            parameters=data.get("metadata", {}).get("parameters", {}),
            output_path=output_path,
            pre_rendered_candidates=candidate_dicts,
        )

    click.secho(f"✅ Report written to: {output}", fg="green")


@cli.command("check-llm")
def check_llm() -> None:
    """Check available LLM backends."""
    click.echo(f"🤖 DogNeo v{__version__} — LLM Backend Check")
    click.echo()

    from dogneo.llm.cli_wrapper import check_cli_availability
    from dogneo.config import LLMConfig

    # CLI backends
    click.echo("📡 CLI Backends (Tier 1 — Free):")
    cli_status = check_cli_availability()
    for tool, available in cli_status.items():
        icon = "✅" if available else "❌"
        click.echo(f"   {icon} {tool}")

    # Local backends
    click.echo("\n💾 Local Backends (Tier 2 — Offline):")
    config = LLMConfig()
    if config.local_model_path:
        exists = Path(config.local_model_path).exists()
        icon = "✅" if exists else "❌"
        click.echo(f"   {icon} {config.local_model_path}")
    else:
        click.echo("   ⚪ No local model configured")

    # Cloud backends
    click.echo("\n☁️  Cloud Backends (Tier 3 — API):")
    import os
    for key, name in [
        ("OPENAI_API_KEY", "OpenAI"),
        ("ANTHROPIC_API_KEY", "Anthropic"),
        ("GOOGLE_API_KEY", "Google Gemini"),
    ]:
        icon = "✅" if os.environ.get(key) else "❌"
        click.echo(f"   {icon} {name} ({key})")


@cli.command("version")
def version_cmd() -> None:
    """Show version and disclaimer."""
    click.echo(f"DogNeo v{__version__}")
    click.echo(f"⚠️  {RUO_DISCLAIMER}")


def main() -> None:
    """CLI entrypoint."""
    cli()


if __name__ == "__main__":
    main()
