"""Report generator for neoantigen analysis results.

Produces HTML and Markdown reports using Jinja2 templates,
with optional LLM-assisted commentary.
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dogneo import RUO_DISCLAIMER, __version__
from dogneo.core.ranking import NeoantigenCandidate
from dogneo.llm.router import LLMRouter, TaskType

logger = logging.getLogger(__name__)


# Inline HTML template (avoids external file dependency for simple usage)
HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>DogNeo Report — {{ sample_id }}</title>
    <style>
        :root { --bg: #0f172a; --surface: #1e293b; --text: #e2e8f0;
                --accent: #38bdf8; --warn: #f97316; --muted: #94a3b8; }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body { font-family: 'Inter', system-ui, sans-serif; background: var(--bg);
               color: var(--text); line-height: 1.6; padding: 2rem; }
        .container { max-width: 1100px; margin: auto; }
        h1 { color: var(--accent); font-size: 1.8rem; margin-bottom: 0.5rem; }
        h2 { color: var(--accent); font-size: 1.3rem; margin: 1.5rem 0 0.5rem;
             border-bottom: 1px solid var(--surface); padding-bottom: 0.3rem; }
        .disclaimer { background: #7c2d12; color: #fed7aa; padding: 1rem;
                       border-radius: 8px; margin: 1rem 0; font-weight: 600; }
        .meta { color: var(--muted); font-size: 0.85rem; margin-bottom: 1rem; }
        table { width: 100%; border-collapse: collapse; margin: 1rem 0;
                background: var(--surface); border-radius: 8px; overflow: hidden; }
        th { background: #334155; color: var(--accent); padding: 0.7rem 0.5rem;
             text-align: left; font-size: 0.8rem; text-transform: uppercase; }
        td { padding: 0.5rem; border-top: 1px solid #334155; font-size: 0.85rem; }
        tr:hover td { background: #334155; }
        .score-high { color: #4ade80; font-weight: 700; }
        .score-med { color: #facc15; }
        .score-low { color: var(--muted); }
        .summary { background: var(--surface); padding: 1.2rem;
                    border-radius: 8px; margin: 1rem 0; }
        .footer { color: var(--muted); font-size: 0.75rem; margin-top: 2rem;
                   text-align: center; }
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 DogNeo — Neoantigen Analysis Report</h1>
    <p class="meta">
        Sample: <strong>{{ sample_id }}</strong> |
        Generated: {{ timestamp }} |
        DogNeo v{{ version }}
    </p>

    <div class="disclaimer">⚠️ {{ ruo_disclaimer }}</div>

    {% if ai_summary %}
    <h2>AI-Assisted Summary</h2>
    <div class="summary">{{ ai_summary }}</div>
    {% endif %}

    <h2>Pipeline Parameters</h2>
    <table>
        <tr><th>Parameter</th><th>Value</th></tr>
        {% for key, value in parameters.items() %}
        <tr><td>{{ key }}</td><td>{{ value }}</td></tr>
        {% endfor %}
    </table>

    <h2>Top Candidates ({{ candidates|length }} total)</h2>
    <table>
        <tr>
            <th>#</th><th>Gene</th><th>Mutation</th><th>Peptide</th>
            <th>Allele</th><th>Affinity (nM)</th><th>TPM</th>
            <th>VAF</th><th>Score</th>
        </tr>
        {% for c in candidates %}
        <tr>
            <td>{{ c.rank }}</td>
            <td>{{ c.gene }}</td>
            <td>{{ c.mutation }}</td>
            <td><code>{{ c.mutant_peptide }}</code></td>
            <td>{{ c.allele }}</td>
            <td>{{ "%.1f"|format(c.binding_affinity_nm) }}</td>
            <td>{{ "%.1f"|format(c.expression_tpm) }}</td>
            <td>{{ "%.2f"|format(c.vaf) }}</td>
            <td class="{{ 'score-high' if c.composite_score > 0.7 else 'score-med' if c.composite_score > 0.4 else 'score-low' }}">
                {{ "%.4f"|format(c.composite_score) }}
            </td>
        </tr>
        {% endfor %}
    </table>

    {% if alleles %}
    <h2>DLA Alleles</h2>
    <ul>
        {% for allele in alleles %}
        <li>{{ allele }}</li>
        {% endfor %}
    </ul>
    {% endif %}

    <div class="footer">
        DogNeo v{{ version }} | {{ ruo_disclaimer }}
    </div>
</div>
</body>
</html>
"""


class ReportGenerator:
    """Generate analysis reports from ranked neoantigen candidates.

    Supports HTML and Markdown output. Optionally uses an LLM router
    for AI-assisted summary and analysis.

    Args:
        llm_router: Optional LLM router for AI commentary.
    """

    def __init__(self, llm_router: LLMRouter | None = None):
        self.llm_router = llm_router

    def generate_html(
        self,
        candidates: list[NeoantigenCandidate],
        sample_id: str,
        parameters: dict[str, Any] | None = None,
        alleles: list[str] | None = None,
        output_path: str | Path | None = None,
        top_n: int = 50,
        pre_rendered_candidates: list[dict] | None = None,
    ) -> str:
        """Generate an HTML report from ranked candidates.

        Args:
            candidates: Ranked neoantigen candidates.
            sample_id: Sample identifier.
            parameters: Pipeline parameters dict.
            alleles: DLA alleles used.
            output_path: Output file path (optional).
            top_n: Number of candidates to show.

        Returns:
            Rendered HTML string.
        """
        from jinja2 import Template

        # Prepare candidate dicts for template
        if pre_rendered_candidates is not None:
            candidate_dicts = pre_rendered_candidates[:top_n]
        else:
            candidate_dicts = [c.to_dict() for c in candidates[:top_n]]

        # Generate AI summary if router available
        ai_summary = ""
        if self.llm_router:
            try:
                from dogneo.llm.prompts import build_rank_analysis_prompt
                prompt = build_rank_analysis_prompt(candidate_dicts, min(10, len(candidate_dicts)))
                ai_summary = self.llm_router.generate(prompt, task_type=TaskType.SUMMARIZE)
            except Exception as e:
                logger.warning("AI summary generation failed: %s", e)

        # Render template
        template = Template(HTML_TEMPLATE)
        html = template.render(
            sample_id=sample_id,
            timestamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            version=__version__,
            ruo_disclaimer=RUO_DISCLAIMER,
            candidates=candidate_dicts,
            parameters=parameters or {},
            alleles=alleles or [],
            ai_summary=ai_summary,
        )

        if output_path:
            Path(output_path).write_text(html, encoding="utf-8")
            logger.info("HTML report written to: %s", output_path)

        return html

    def generate_markdown(
        self,
        candidates: list[NeoantigenCandidate],
        sample_id: str,
        top_n: int = 20,
        parameters: dict[str, Any] | None = None,
        output_path: str | Path | None = None,
        pre_rendered_candidates: list[dict] | None = None,
    ) -> str:
        """Generate a Markdown summary of top candidates.

        Args:
            candidates: Ranked neoantigen candidates.
            sample_id: Sample identifier.
            top_n: Number of candidates to include.

        Returns:
            Markdown text.
        """
        lines = [
            f"# Neoantigen Analysis — {sample_id}",
            "",
            f"> {RUO_DISCLAIMER}",
            "",
            f"**Generated**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}  ",
            f"**DogNeo**: v{__version__}",
            "",
            "## Top Candidates",
            "",
            "| # | Gene | Mutation | Peptide | Allele | Affinity | TPM | Score |",
            "|---|------|----------|---------|--------|----------|-----|-------|",
        ]

        if pre_rendered_candidates is not None:
            candidate_dicts = pre_rendered_candidates[:top_n]
        else:
            candidate_dicts = [c.to_dict() for c in candidates[:top_n]]

        for d in candidate_dicts:
            lines.append(
                f"| {d.get('rank', '-')} | {d.get('gene', '')} | {d.get('mutation', '')} | "
                f"`{d.get('mutant_peptide', '')}` | {d.get('allele', '')} | "
                f"{float(d.get('binding_affinity_nm', 0)):.1f} | {float(d.get('expression_tpm', 0)):.1f} | "
                f"{float(d.get('composite_score', 0)):.4f} |"
            )

        total = len(pre_rendered_candidates) if pre_rendered_candidates else len(candidates)
        lines.extend(["", f"*Total candidates: {total}*"])
        md = "\n".join(lines)

        if output_path:
            Path(output_path).write_text(md, encoding="utf-8")
            logger.info("Markdown report written to: %s", output_path)

        return md
