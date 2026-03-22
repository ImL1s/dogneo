"""Tests for ReferenceDataManager — canine reference data download & caching.

TDD RED phase: these tests define the expected API surface.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


class TestReferenceDataManager:
    """Unit tests for ReferenceDataManager (offline, mocked HTTP)."""

    def test_import(self):
        """Manager can be imported."""
        from dogneo.data.manager import ReferenceDataManager
        assert ReferenceDataManager is not None

    def test_default_cache_dir(self):
        """Default cache dir is ~/.dogneo/data/."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager()
        assert mgr.cache_dir == Path.home() / ".dogneo" / "data"

    def test_custom_cache_dir(self, tmp_path):
        """Cache dir can be overridden."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager(cache_dir=tmp_path / "custom")
        assert mgr.cache_dir == tmp_path / "custom"

    def test_get_proteome_path_not_setup(self, tmp_path):
        """Returns None when proteome hasn't been downloaded yet."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager(cache_dir=tmp_path / "empty")
        assert mgr.get_proteome_path() is None

    def test_get_proteome_path_after_setup(self, tmp_path):
        """Returns the FASTA path when proteome exists in cache."""
        from dogneo.data.manager import ReferenceDataManager

        cache = tmp_path / "cache"
        cache.mkdir()
        # Simulate a cached proteome
        fa = cache / "CanFam3.1.pep.all.fa"
        fa.write_text(">ENSCAFP00000000001\nMKTAYIAKQ\n")

        mgr = ReferenceDataManager(cache_dir=cache)
        result = mgr.get_proteome_path()
        assert result == fa
        assert result.exists()

    def test_get_dla_alleles_path_bundled(self, tmp_path):
        """DLA alleles path points to bundled package data."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager(cache_dir=tmp_path)
        path = mgr.get_dla_alleles_path()
        assert path.exists()
        assert path.name == "dla_alleles.txt"

    def test_status_not_setup(self, tmp_path):
        """Status reports proteome as not downloaded."""
        from dogneo.data.manager import ReferenceDataManager

        mgr = ReferenceDataManager(cache_dir=tmp_path / "nope")
        status = mgr.status()
        assert status["proteome_ready"] is False
        assert status["dla_alleles_ready"] is True  # bundled

    def test_status_after_setup(self, tmp_path):
        """Status reports proteome as ready after setup."""
        from dogneo.data.manager import ReferenceDataManager

        cache = tmp_path / "cache"
        cache.mkdir()
        fa = cache / "CanFam3.1.pep.all.fa"
        fa.write_text(">SEQ1\nMKTAYIAKQ\n")

        mgr = ReferenceDataManager(cache_dir=cache)
        status = mgr.status()
        assert status["proteome_ready"] is True
        assert status["proteome_proteins"] > 0

    def test_setup_downloads_and_decompresses(self, tmp_path):
        """setup() downloads gzipped FASTA and decompresses it."""
        from dogneo.data.manager import ReferenceDataManager

        # Create fake .fa.gz content
        fasta_content = b">ENSCAFP00000000001\nMKTAYIAKQ\n>ENSCAFP00000000002\nMWKKL\n"
        gz_content = gzip.compress(fasta_content)

        # Mock the HTTP response
        mock_response = MagicMock()
        mock_response.ok = True
        mock_response.headers = {"content-length": str(len(gz_content))}
        mock_response.iter_content = lambda chunk_size: [gz_content]
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = MagicMock(return_value=False)

        cache = tmp_path / "dl_cache"

        with patch("dogneo.data.manager.requests.get", return_value=mock_response):
            mgr = ReferenceDataManager(cache_dir=cache)
            path = mgr.setup()

        assert path.exists()
        assert path.suffix == ".fa"
        content = path.read_text()
        assert ">ENSCAFP00000000001" in content
        assert "MKTAYIAKQ" in content

    def test_setup_skips_if_cached(self, tmp_path):
        """setup() does not re-download if cached file exists."""
        from dogneo.data.manager import ReferenceDataManager

        cache = tmp_path / "cache"
        cache.mkdir()
        fa = cache / "CanFam3.1.pep.all.fa"
        fa.write_text(">SEQ1\nMKTAYIAKQ\n")

        with patch("dogneo.data.manager.requests.get") as mock_get:
            mgr = ReferenceDataManager(cache_dir=cache)
            path = mgr.setup()
            mock_get.assert_not_called()

        assert path == fa

    def test_setup_force_redownloads(self, tmp_path):
        """setup(force=True) re-downloads even if cached."""
        from dogneo.data.manager import ReferenceDataManager

        cache = tmp_path / "cache"
        cache.mkdir()
        fa = cache / "CanFam3.1.pep.all.fa"
        fa.write_text(">OLD\nMKK\n")

        new_content = b">NEW\nMKTAYIAKQ\n"
        gz_content = gzip.compress(new_content)

        mock_response = MagicMock()
        mock_response.ok = True
        mock_response.headers = {"content-length": str(len(gz_content))}
        mock_response.iter_content = lambda chunk_size: [gz_content]
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = MagicMock(return_value=False)

        with patch("dogneo.data.manager.requests.get", return_value=mock_response):
            mgr = ReferenceDataManager(cache_dir=cache)
            path = mgr.setup(force=True)

        assert ">NEW" in path.read_text()

    def test_setup_http_failure_raises(self, tmp_path):
        """setup() raises on HTTP failure."""
        from dogneo.data.manager import ReferenceDataManager

        mock_response = MagicMock()
        mock_response.ok = False
        mock_response.status_code = 404
        mock_response.text = "Not Found"
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = MagicMock(return_value=False)

        with patch("dogneo.data.manager.requests.get", return_value=mock_response):
            mgr = ReferenceDataManager(cache_dir=tmp_path / "fail")
            with pytest.raises(RuntimeError, match="download"):
                mgr.setup()

    def test_count_proteins(self, tmp_path):
        """_count_proteins correctly counts FASTA entries."""
        from dogneo.data.manager import ReferenceDataManager

        fa = tmp_path / "test.fa"
        fa.write_text(">P1\nMKK\n>P2\nMWW\n>P3\nMAA\n")

        mgr = ReferenceDataManager(cache_dir=tmp_path)
        assert mgr._count_proteins(fa) == 3

    def test_env_override_cache_dir(self, tmp_path, monkeypatch):
        """DOGNEO_DATA_DIR env var overrides default cache dir."""
        from dogneo.data.manager import ReferenceDataManager

        monkeypatch.setenv("DOGNEO_DATA_DIR", str(tmp_path / "env_dir"))
        mgr = ReferenceDataManager()
        assert mgr.cache_dir == tmp_path / "env_dir"
