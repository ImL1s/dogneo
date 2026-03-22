"""Tests for IEDB API client with caching and rate limiting.

TDD Phase 2: IEDB binding prediction for out-of-box scored results.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from dogneo.core.binding import BindingPrediction


class TestIEDBClient:
    """Test IEDBClient core functionality."""

    def test_client_can_be_instantiated(self):
        from dogneo.core.iedb_client import IEDBClient

        client = IEDBClient()
        assert client is not None

    def test_client_accepts_custom_cache_dir(self, tmp_path):
        from dogneo.core.iedb_client import IEDBClient

        client = IEDBClient(cache_dir=tmp_path / "my_cache")
        assert client.cache_dir == tmp_path / "my_cache"

    def test_predict_batch_returns_binding_predictions(self, tmp_path):
        """Mock the HTTP call, verify we get BindingPrediction objects back."""
        from dogneo.core.iedb_client import IEDBClient

        mock_response_text = (
            "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tpercentile_rank\tic50\n"
            "DLA-88*001:01\t1\t1\t9\t9\tRAIVGAPPS\tnetmhcpan_el\t1.5\t125.0\n"
            "DLA-88*001:01\t2\t1\t9\t9\tYLVSGAPPS\tnetmhcpan_el\t0.3\t42.0\n"
        )

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = mock_response_text

        client = IEDBClient(cache_dir=tmp_path / "cache")

        with patch("requests.post", return_value=mock_resp):
            results = client.predict_batch(
                peptides=["RAIVGAPPS", "YLVSGAPPS"],
                alleles=["DLA-88*001:01"],
            )

        assert len(results) == 2
        assert all(isinstance(r, BindingPrediction) for r in results)
        assert results[0].peptide_sequence == "RAIVGAPPS"
        assert results[0].allele == "DLA-88*001:01"
        assert results[0].affinity_nm == 125.0
        assert results[0].percentile_rank == 1.5
        assert results[0].tool == "iedb"
        assert results[1].affinity_nm == 42.0

    def test_predict_batch_caches_results(self, tmp_path):
        """After a successful call, results should be cached to disk."""
        from dogneo.core.iedb_client import IEDBClient

        mock_response_text = (
            "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tpercentile_rank\tic50\n"
            "DLA-88*001:01\t1\t1\t9\t9\tRAIVGAPPS\tnetmhcpan_el\t1.5\t125.0\n"
        )

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = mock_response_text

        cache_dir = tmp_path / "cache"
        client = IEDBClient(cache_dir=cache_dir)

        with patch("requests.post", return_value=mock_resp):
            client.predict_batch(
                peptides=["RAIVGAPPS"],
                alleles=["DLA-88*001:01"],
            )

        # Cache directory should exist and have files
        assert cache_dir.exists()
        cache_files = list(cache_dir.glob("*.json"))
        assert len(cache_files) > 0

    def test_predict_batch_uses_cache_on_second_call(self, tmp_path):
        """Second call with same inputs should use cache, not HTTP."""
        from dogneo.core.iedb_client import IEDBClient

        mock_response_text = (
            "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tpercentile_rank\tic50\n"
            "DLA-88*001:01\t1\t1\t9\t9\tRAIVGAPPS\tnetmhcpan_el\t1.5\t125.0\n"
        )

        mock_resp = MagicMock()
        mock_resp.ok = True
        mock_resp.text = mock_response_text

        client = IEDBClient(cache_dir=tmp_path / "cache")

        with patch("requests.post", return_value=mock_resp) as mock_post:
            # First call — hits API
            result1 = client.predict_batch(["RAIVGAPPS"], ["DLA-88*001:01"])
            assert mock_post.call_count == 1

            # Second call — should use cache
            result2 = client.predict_batch(["RAIVGAPPS"], ["DLA-88*001:01"])
            assert mock_post.call_count == 1  # Still 1, no new call

        assert len(result1) == len(result2)
        assert result1[0].affinity_nm == result2[0].affinity_nm

    def test_predict_batch_handles_offline(self, tmp_path):
        """When network is unavailable, should raise ConnectionError."""
        from dogneo.core.iedb_client import IEDBClient

        client = IEDBClient(cache_dir=tmp_path / "cache")

        import requests
        with patch("requests.post", side_effect=requests.ConnectionError("offline")):
            with pytest.raises(requests.ConnectionError):
                client.predict_batch(["RAIVGAPPS"], ["DLA-88*001:01"])

    def test_predict_batch_empty_inputs(self, tmp_path):
        """Empty peptides or alleles should return empty list."""
        from dogneo.core.iedb_client import IEDBClient

        client = IEDBClient(cache_dir=tmp_path / "cache")
        assert client.predict_batch([], ["DLA-88*001:01"]) == []
        assert client.predict_batch(["RAIVGAPPS"], []) == []


class TestBindingFallback:
    """Test auto-fallback chain: netmhcpan → iedb → none."""

    def test_resolve_binding_tool_returns_estimator_when_no_netmhcpan(self):
        from dogneo.app.rank_pipeline import _resolve_binding_tool

        with patch("shutil.which", return_value=None):
            tool = _resolve_binding_tool("auto")
        assert tool == "estimator"

    def test_resolve_binding_tool_returns_netmhcpan_when_available(self):
        from dogneo.app.rank_pipeline import _resolve_binding_tool

        with patch("shutil.which", return_value="/usr/local/bin/netMHCpan"):
            tool = _resolve_binding_tool("auto")
        assert tool == "netmhcpan"

    def test_resolve_binding_tool_returns_explicit_choice(self):
        from dogneo.app.rank_pipeline import _resolve_binding_tool

        assert _resolve_binding_tool("iedb") == "iedb"
        assert _resolve_binding_tool("none") == "none"
        assert _resolve_binding_tool("netmhcpan") == "netmhcpan"
