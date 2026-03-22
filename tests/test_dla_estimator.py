"""Tests for built-in DLA binding estimator.

Provides out-of-box scoring when no external predictor is available.
Based on published DLA-88 binding motifs.
"""
from __future__ import annotations

import pytest


class TestDLAEstimator:

    def test_estimator_returns_binding_predictions(self):
        from dogneo.core.dla_estimator import estimate_binding

        results = estimate_binding(
            peptides=["RAIVGAPPS", "YLVSGAPPS"],
            alleles=["DLA-88*001:01"],
        )
        assert len(results) == 2
        assert all(r.tool == "dogneo-estimator" for r in results)

    def test_estimated_affinity_is_numeric(self):
        from dogneo.core.dla_estimator import estimate_binding

        results = estimate_binding(["RAIVGAPPS"], ["DLA-88*001:01"])
        assert results[0].affinity_nm > 0
        assert results[0].percentile_rank >= 0

    def test_anchor_residues_improve_score(self):
        """Peptides with preferred anchor residues should score better."""
        from dogneo.core.dla_estimator import estimate_binding

        # DLA-88 prefers hydrophobic anchors at P2 and P9
        good_peptide = "KLVFFAEDV"  # L at P2, V at P9 — good anchors
        bad_peptide = "KEEEEEEEE"   # E (charged) at anchors — poor

        good = estimate_binding([good_peptide], ["DLA-88*001:01"])[0]
        bad = estimate_binding([bad_peptide], ["DLA-88*001:01"])[0]

        # Lower affinity = better binding
        assert good.affinity_nm < bad.affinity_nm

    def test_empty_inputs(self):
        from dogneo.core.dla_estimator import estimate_binding

        assert estimate_binding([], ["DLA-88*001:01"]) == []
        assert estimate_binding(["RAIVGAPPS"], []) == []

    def test_multiple_alleles(self):
        from dogneo.core.dla_estimator import estimate_binding

        results = estimate_binding(
            ["RAIVGAPPS"],
            ["DLA-88*001:01", "DLA-88*002:01"],
        )
        assert len(results) == 2

    def test_different_peptide_lengths(self):
        from dogneo.core.dla_estimator import estimate_binding

        for length in [8, 9, 10, 11]:
            pep = "A" * length
            results = estimate_binding([pep], ["DLA-88*001:01"])
            assert len(results) == 1
            assert results[0].affinity_nm > 0
