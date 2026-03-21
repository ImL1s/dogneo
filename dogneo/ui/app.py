"""DogNeo Streamlit UI — Interactive Neoantigen Analysis Dashboard.

Launch with: dogneo ui  (or: streamlit run dogneo/ui/app.py)
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import streamlit as st

from dogneo import __version__, RUO_DISCLAIMER

st.set_page_config(
    page_title=f"DogNeo v{__version__}",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def _run_demo():
    """Run demo pipeline and cache results."""
    from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline
    from dogneo.data.manager import ReferenceDataManager
    import dogneo.data as _data_pkg

    demo_dir = Path(_data_pkg.__file__).parent / "demo"
    vcf_path = demo_dir / "canine_osteosarcoma.vcf"
    expr_path = demo_dir / "expression.sf"

    mgr = ReferenceDataManager()
    proteome = mgr.get_proteome_path()
    if not proteome:
        st.error("Proteome not found. Run `dogneo setup` first.")
        return None

    inp = RankInput(
        vcf_path=vcf_path,
        expression_path=expr_path,
        sample_id="CANINE_OSA_DEMO",
        protein_db_path=proteome,
        binding_tool="none",
        formats=[],
    )
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        return run_rank_pipeline(inp, Path(tmpdir))


def _run_uploaded(vcf_file, expr_file, alleles_str, binding_tool):
    """Run pipeline on uploaded files."""
    from dogneo.app.rank_pipeline import RankInput, run_rank_pipeline

    import tempfile
    tmpdir = tempfile.mkdtemp()  # Managed by Streamlit session lifecycle
    vcf_path = Path(tmpdir) / "uploaded.vcf"
    vcf_file.seek(0)
    vcf_path.write_bytes(vcf_file.read())

    expr_path = None
    if expr_file:
        expr_path = Path(tmpdir) / "expression.sf"
        expr_file.seek(0)
        expr_path.write_bytes(expr_file.read())

    alleles = [a.strip() for a in alleles_str.split(",") if a.strip()] if alleles_str else []

    inp = RankInput(
        vcf_path=vcf_path,
        expression_path=expr_path,
        sample_id="UPLOADED",
        alleles=alleles,
        binding_tool=binding_tool,
        formats=[],
    )
    return run_rank_pipeline(inp, Path(tmpdir))


def main():
    # Sidebar
    st.sidebar.title("🧬 DogNeo")
    st.sidebar.caption(f"v{__version__}")
    st.sidebar.warning(RUO_DISCLAIMER)

    page = st.sidebar.radio(
        "Navigate",
        ["Upload / Demo", "Ranking", "Visualizations", "Report"],
    )

    # ---- Page: Upload / Demo ----
    if page == "Upload / Demo":
        st.title("🧬 DogNeo — Canine Neoantigen Analysis")
        st.markdown("From somatic variants to ranked vaccine candidates — purpose-built for dogs.")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Run Demo")
            st.markdown("Analyze bundled canine osteosarcoma data (8 mutations, 6 coding).")
            if st.button("Run Demo", type="primary", use_container_width=True):
                with st.spinner("Running demo pipeline..."):
                    result = _run_demo()
                if result:
                    st.session_state["result"] = result
                    st.success(f"Demo complete! {len(result.candidates)} candidates.")

        with col2:
            st.subheader("Upload Your Data")
            vcf_file = st.file_uploader("VCF file (annotated)", type=["vcf"])
            expr_file = st.file_uploader("Expression file (optional)", type=["sf", "tsv"])
            alleles_str = st.text_input("DLA alleles (comma-separated)", placeholder="DLA-88*001:01")
            binding_tool = st.selectbox("Binding predictor", ["auto", "iedb", "none"])

            if st.button("Analyze", use_container_width=True) and vcf_file:
                with st.spinner("Running pipeline..."):
                    result = _run_uploaded(vcf_file, expr_file, alleles_str, binding_tool)
                st.session_state["result"] = result
                st.success(f"Analysis complete! {len(result.candidates)} candidates.")

        # Show summary if result exists
        if "result" in st.session_state:
            result = st.session_state["result"]
            st.divider()
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Variants", result.variants_total)
            c2.metric("Coding", result.variants_coding)
            c3.metric("Peptides", result.peptides_total)
            c4.metric("Candidates", len(result.candidates))

    # ---- Page: Ranking ----
    elif page == "Ranking":
        st.title("Ranking Table")

        if "result" not in st.session_state:
            st.info("Run a demo or upload data first.")
            return

        result = st.session_state["result"]
        candidates = result.candidates

        if not candidates:
            st.warning("No candidates to display.")
            return

        import pandas as pd
        rows = [c.to_dict() for c in candidates]
        df = pd.DataFrame(rows)

        display_cols = [
            "rank", "gene", "mutation", "mutant_peptide", "allele",
            "binding_affinity_nm", "expression_tpm", "vaf", "composite_score",
        ]
        existing_cols = [c for c in display_cols if c in df.columns]
        st.dataframe(
            df[existing_cols],
            use_container_width=True,
            hide_index=True,
        )

        # Export buttons
        col1, col2, col3 = st.columns(3)
        with col1:
            csv_data = df[existing_cols].to_csv(index=False, sep="\t")
            st.download_button("Download TSV", csv_data, "candidates.tsv", "text/tab-separated-values")
        with col2:
            json_data = json.dumps({"candidates": rows}, indent=2)
            st.download_button("Download JSON", json_data, "candidates.json", "application/json")
        with col3:
            fasta_lines = []
            for c in candidates:
                d = c.to_dict()
                fasta_lines.append(f">{d['gene']}_{d['mutation']} rank={d['rank']} score={d['composite_score']:.4f}")
                fasta_lines.append(d["mutant_peptide"])
            st.download_button("Download FASTA", "\n".join(fasta_lines), "candidates.fasta", "text/plain")

    # ---- Page: Visualizations ----
    elif page == "Visualizations":
        st.title("Visualizations")

        if "result" not in st.session_state:
            st.info("Run a demo or upload data first.")
            return

        result = st.session_state["result"]
        candidates = result.candidates

        if not candidates:
            st.warning("No candidates to visualize.")
            return

        from dogneo.ui.charts import (
            score_distribution_chart,
            binding_heatmap,
            score_radar_chart,
            candidate_scatter,
        )

        tab1, tab2, tab3, tab4 = st.tabs(["Score Distribution", "Binding Heatmap", "Radar", "Scatter"])

        with tab1:
            fig = score_distribution_chart(candidates[:30])
            st.plotly_chart(fig, use_container_width=True)

        with tab2:
            fig = binding_heatmap(candidates[:20])
            st.plotly_chart(fig, use_container_width=True)

        with tab3:
            idx = st.selectbox(
                "Select candidate",
                range(min(10, len(candidates))),
                format_func=lambda i: f"#{i+1} {candidates[i].variant.gene} {candidates[i].peptide.mutation}",
            )
            if candidates[idx].score_components:
                fig = score_radar_chart(candidates[idx])
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No score components available (unscored candidate).")

        with tab4:
            col1, col2 = st.columns(2)
            x_axis = col1.selectbox("X axis", ["binding_affinity_nm", "expression_tpm", "vaf", "composite_score"])
            y_axis = col2.selectbox("Y axis", ["expression_tpm", "binding_affinity_nm", "vaf", "composite_score"], index=1)
            fig = candidate_scatter(candidates, x=x_axis, y=y_axis)
            st.plotly_chart(fig, use_container_width=True)

    # ---- Page: Report ----
    elif page == "Report":
        st.title("Analysis Report")

        if "result" not in st.session_state:
            st.info("Run a demo or upload data first.")
            return

        result = st.session_state["result"]

        from dogneo.report.generator import ReportGenerator
        gen = ReportGenerator()
        html = gen.generate_html(
            result.candidates,
            sample_id="DogNeo Analysis",
            parameters={
                "Variants": result.variants_total,
                "Coding": result.variants_coding,
                "Peptides": result.peptides_total,
                "Binding Tool": result.binding_tool_used,
            },
            alleles=result.alleles_used,
        )
        st.components.v1.html(html, height=800, scrolling=True)

        st.download_button("Download HTML Report", html, "report.html", "text/html")


if __name__ == "__main__":
    main()
