"""Prompt templates for neoantigen analysis tasks.

Provides structured prompts for LLM-assisted analysis,
with bilingual (English + Chinese) support.
"""

from __future__ import annotations

from string import Template

# ---------------------------------------------------------------------------
# Prompt Templates
# ---------------------------------------------------------------------------

RANK_ANALYSIS = Template("""You are a tumor immunology expert specializing in canine oncology.

Below is a ranked list of candidate neoantigens from a canine tumor sample:

$candidates_table

For each top candidate, evaluate:
1. Binding affinity strength and clinical relevance
2. Expression level adequacy for immune recognition
3. Self-similarity to wild-type (immunogenicity potential)
4. Overall recommendation (Priority: HIGH / MEDIUM / LOW)

Provide your analysis in a concise, actionable format.
Focus on the top $top_n candidates.

IMPORTANT: This is for RESEARCH USE ONLY. Do not provide clinical or veterinary advice.""")


CANDIDATE_REPORT = Template("""Given a list of top candidate neoantigen peptides with their scores:

$candidates_json

Generate a concise Markdown report with:
1. Executive summary (3-5 sentences)
2. Top candidates table with key metrics
3. Methodology notes
4. Limitations and caveats

Include this disclaimer at the end:
"FOR RESEARCH USE ONLY — NOT FOR CLINICAL OR VETERINARY DIAGNOSTIC USE."

Output in $language.""")


WORKFLOW_SUMMARY = Template("""Summarize the following bioinformatics workflow execution:

Sample ID: $sample_id
Reference Genome: $reference_genome
Steps completed:
$steps_summary

Results:
- Total variants detected: $total_variants
- Coding variants: $coding_variants
- Candidate neoantigens: $total_candidates
- Top candidate: $top_candidate

Generate a concise summary suitable for a lab notebook entry.
Output in $language.""")


LAB_NOTEBOOK = Template("""以下是犬隻腫瘤新抗原分析的流程記錄：

樣本 ID: $sample_id
參考基因組: $reference_genome
分析日期: $analysis_date

流程步驟：
$steps_summary

結果摘要：
- 偵測到的體細胞突變數: $total_variants
- 編碼區突變數: $coding_variants
- 候選新抗原數: $total_candidates
- 排名第一的候選: $top_candidate

請以繁體中文撰寫實驗室筆記格式的摘要，包含：
1. 目的
2. 方法概述
3. 結果
4. 討論與下一步

請注意：此為研究用途，非臨床使用。""")


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

def format_candidates_table(candidates: list[dict], top_n: int = 10) -> str:
    """Format candidate list as a text table for prompt injection.

    Args:
        candidates: List of candidate dicts (from NeoantigenCandidate.to_dict()).
        top_n: Number of top candidates to include.

    Returns:
        Formatted text table.
    """
    header = (
        f"{'Rank':<5} {'Gene':<10} {'Mutation':<15} {'Peptide':<15} "
        f"{'Allele':<20} {'Affinity(nM)':<14} {'TPM':<10} {'Score':<8}"
    )
    lines = [header, "-" * len(header)]

    for c in candidates[:top_n]:
        lines.append(
            f"{c.get('rank', ''):>4}  {c.get('gene', ''):<10} "
            f"{c.get('mutation', ''):<15} {c.get('mutant_peptide', ''):<15} "
            f"{c.get('allele', ''):<20} {c.get('binding_affinity_nm', 0):>12.1f} "
            f"{c.get('expression_tpm', 0):>8.1f} {c.get('composite_score', 0):>7.4f}"
        )

    return "\n".join(lines)


def build_rank_analysis_prompt(
    candidates: list[dict],
    top_n: int = 10,
) -> str:
    """Build a rank analysis prompt from candidate data.

    Args:
        candidates: Candidate dicts.
        top_n: Number to analyze.

    Returns:
        Formatted prompt string.
    """
    table = format_candidates_table(candidates, top_n)
    return RANK_ANALYSIS.substitute(
        candidates_table=table,
        top_n=top_n,
    )


def build_report_prompt(
    candidates_json: str,
    language: str = "English",
) -> str:
    """Build a report generation prompt.

    Args:
        candidates_json: JSON string of candidates.
        language: Output language.

    Returns:
        Formatted prompt string.
    """
    return CANDIDATE_REPORT.substitute(
        candidates_json=candidates_json,
        language=language,
    )


def build_workflow_summary_prompt(
    sample_id: str,
    reference_genome: str,
    steps_summary: str,
    total_variants: int,
    coding_variants: int,
    total_candidates: int,
    top_candidate: str,
    language: str = "English",
) -> str:
    """Build a workflow summary prompt.

    Args:
        sample_id: Sample identifier.
        reference_genome: Reference genome used.
        steps_summary: Text summary of completed steps.
        total_variants: Number of variants found.
        coding_variants: Number of coding variants.
        total_candidates: Number of neoantigen candidates.
        top_candidate: Description of top candidate.
        language: Output language.

    Returns:
        Formatted prompt string.
    """
    return WORKFLOW_SUMMARY.substitute(
        sample_id=sample_id,
        reference_genome=reference_genome,
        steps_summary=steps_summary,
        total_variants=total_variants,
        coding_variants=coding_variants,
        total_candidates=total_candidates,
        top_candidate=top_candidate,
        language=language,
    )
