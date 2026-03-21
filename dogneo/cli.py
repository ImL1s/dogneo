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
@click.option("--protein-db", type=click.Path(exists=True),
              help="Canine protein FASTA database (e.g., CanFam3.1 proteome).")
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
    protein_db: str | None,
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

    # Parse alleles — auto-load bundled DLA alleles if none provided
    allele_list = [a.strip() for a in alleles.split(",") if a.strip()] if alleles else []
    if not allele_list:
        from dogneo.data.manager import ReferenceDataManager
        mgr = ReferenceDataManager()
        dla_path = mgr.get_dla_alleles_path()
        if dla_path.exists():
            allele_list = [l.strip() for l in open(dla_path)
                          if l.strip() and not l.strip().startswith('#')][:4]
            click.echo(f"   Auto-loaded {len(allele_list)} DLA alleles from bundled data")

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

    pdb = ProteinDatabase()

    # Auto-detect proteome: explicit --protein-db > cached > none
    effective_protein_db = protein_db
    if not effective_protein_db:
        from dogneo.data.manager import ReferenceDataManager
        mgr = ReferenceDataManager()
        cached = mgr.get_proteome_path()
        if cached:
            effective_protein_db = str(cached)
            click.echo(f"🧬 Auto-detected cached proteome: {cached.name}")

    if effective_protein_db:
        click.echo(f"🧬 Loading protein database: {Path(effective_protein_db).name}")
        pdb.load_fasta(effective_protein_db)
    peptides_by_variant: dict[str, list] = {}
    predictions_by_peptide: dict[str, list] = {}

    for v in coding:
        peps = generate_peptides(v, pdb, lengths=mhci_lens)
        if peps:
            peptides_by_variant[v.variant_id] = peps

    total_peptides = sum(len(peps) for peps in peptides_by_variant.values())
    click.echo(f"   {total_peptides} peptides from {len(peptides_by_variant)} variants")

    if not peptides_by_variant:
        click.secho(
            "⚠️  No peptides generated. Try:\n"
            "   • Run `dogneo setup` to download the CanFam3.1 proteome\n"
            "   • Or provide --protein-db /path/to/proteome.fa",
            fg="yellow",
        )

    # Step 4: Ranking
    click.echo("📊 Scoring and ranking candidates...")

    # If no binding predictions are available (standalone mode),
    # create candidates with placeholder binding for export.
    if not predictions_by_peptide and peptides_by_variant:
        click.echo("   ℹ️  No binding predictor — creating unscored candidates")
        click.echo("   💡 Run NetMHCpan on the exported FASTA, then re-rank.")
        candidates_list: list = []
        for variant in coding:
            peptides = peptides_by_variant.get(variant.variant_id, [])
            for peptide in peptides:
                # Create placeholder binding for each allele
                for allele in (allele_list or ["DLA-88*unknown"]):
                    candidates_list.append(NeoantigenCandidate(
                        variant=variant,
                        peptide=peptide,
                        binding=BindingPrediction(
                            peptide_sequence=peptide.mut_sequence,
                            allele=allele,
                            affinity_nm=float('nan'),  # unknown — run NetMHCpan
                            percentile_rank=float('nan'),
                            tool="pending",
                            mhc_class=1,
                        ),
                        expression_tpm=variant.expression_tpm,
                    ))
        ranked = candidates_list
        click.echo(f"   {len(ranked)} unscored candidates")
    else:
        candidates_list = build_candidates(coding, peptides_by_variant, predictions_by_peptide)
        if allele_list and candidates_list:
            ranked = rank_candidates(candidates_list)
            click.echo(f"   {len(ranked)} candidates ranked")
        else:
            ranked = candidates_list
            click.echo(f"   {len(ranked)} candidates (unranked — no alleles or binding data)")

    # Step 5: Export
    if "tsv" in format_list:
        export_tsv(ranked, outdir / "candidates.tsv")
    if "json" in format_list:
        export_json(ranked, outdir / "candidates.json", sample_id)
    if "fasta" in format_list:
        export_fasta(ranked, outdir / "candidates.fasta")
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
            ranked, sample_id,
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


@cli.command()
@click.option("--force", is_flag=True, help="Re-download even if cached.")
def setup(force: bool) -> None:
    """Download reference data (CanFam3.1 proteome, DLA alleles).

    Downloads the canine proteome from Ensembl FTP and caches it locally.
    Only needs to be run once. Use --force to re-download.
    """
    click.echo(f"🧬 DogNeo v{__version__} — Setting up reference data")

    from dogneo.data.manager import ReferenceDataManager

    mgr = ReferenceDataManager()
    proteome = mgr.get_proteome_path()

    if proteome and not force:
        n = mgr._count_proteins(proteome)
        click.echo(f"✅ Proteome already cached: {proteome} ({n:,} proteins)")
    else:
        click.echo("📥 Downloading CanFam3.1 proteome from Ensembl...")
        path = mgr.setup(force=force)
        n = mgr._count_proteins(path)
        click.secho(f"✅ Proteome saved: {path} ({n:,} proteins)", fg="green")

    # DLA alleles (bundled)
    dla_path = mgr.get_dla_alleles_path()
    n_alleles = sum(1 for line in open(dla_path) if line.strip())
    click.echo(f"✅ DLA alleles: bundled ({n_alleles} alleles from IPD-MHC)")

    click.secho("\n✅ Setup complete! Run `dogneo demo` to test the pipeline.", fg="green")


@cli.command()
@click.option("--output-dir", "-o", default="dogneo_demo_results",
              help="Output directory for demo results.")
@click.option("--binding", type=click.Choice(["none", "iedb"]),
              default="none", help="Binding predictor (iedb = free online API).")
def demo(output_dir: str, binding: str) -> None:
    """Run a demo pipeline with bundled canine osteosarcoma data.

    Uses built-in demo VCF, expression, and DLA allele data.
    Proteome must be downloaded first via `dogneo setup`.
    """
    click.echo(f"⚠️  {RUO_DISCLAIMER}")
    click.echo(f"🧬 DogNeo v{__version__} — Running demo pipeline")

    from dogneo.data.manager import ReferenceDataManager

    mgr = ReferenceDataManager()
    proteome = mgr.get_proteome_path()

    if not proteome:
        click.secho(
            "❌ Proteome not found. Run `dogneo setup` first to download "
            "the CanFam3.1 reference data.",
            fg="red",
        )
        sys.exit(1)

    # Locate bundled demo data (inside the package: dogneo/data/demo/)
    import dogneo.data as _data_pkg
    demo_dir = Path(_data_pkg.__file__).parent / "demo"
    vcf_path = demo_dir / "canine_osteosarcoma.vcf"
    expr_path = demo_dir / "expression.sf"
    alleles_path = mgr.get_dla_alleles_path()

    if not vcf_path.exists():
        click.secho("❌ Demo VCF not found. Package may be incomplete.", fg="red")
        sys.exit(1)

    click.echo(f"📂 Using bundled demo data: {demo_dir}")
    click.echo(f"📂 Using reference: {proteome}")

    # Load alleles
    alleles = [line.strip() for line in open(alleles_path)
               if line.strip() and not line.strip().startswith('#')]
    allele_str = ",".join(alleles[:4])  # Use first 4 for demo speed

    # Invoke the rank command programmatically
    from click.testing import CliRunner
    runner = CliRunner()
    rank_args = [
        "rank",
        "--vcf", str(vcf_path),
        "--expression", str(expr_path),
        "--sample-id", "CANINE_OSA_DEMO",
        "--output-dir", output_dir,
        "--alleles", allele_str,
        "--protein-db", str(proteome),
        "--formats", "tsv,json,fasta",
    ]

    result = runner.invoke(cli, rank_args, catch_exceptions=False)
    click.echo(result.output)

    if result.exit_code != 0:
        click.secho("❌ Demo pipeline failed.", fg="red")
        sys.exit(result.exit_code)

    out_path = Path(output_dir) / "CANINE_OSA_DEMO"
    click.secho(f"\n✅ Demo complete! Results: {out_path}", fg="green")
    click.echo("   📄 candidates.tsv  — Tab-separated candidates")
    click.echo("   📄 candidates.json — Structured JSON with metadata")
    click.echo("   📄 candidates.fasta — Peptide sequences for wet lab")


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
