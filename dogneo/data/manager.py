"""Reference data manager for DogNeo.

Handles downloading, caching, and validation of canine reference data
(CanFam3.1 proteome, DLA alleles) from public databases.

Inspired by PyEnsembl's `pyensembl install` and pVACtools' `download_example_data`.
"""
from __future__ import annotations

import gzip
import logging
import os
import shutil
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# Ensembl FTP URL for CanFam3.1 proteome (release 104)
PROTEOME_URL = (
    "https://ftp.ensembl.org/pub/release-104/fasta/"
    "canis_lupus_familiaris/pep/Canis_lupus_familiaris.CanFam3.1.pep.all.fa.gz"
)

PROTEOME_FILENAME = "CanFam3.1.pep.all.fa"
PROTEOME_GZ_FILENAME = "CanFam3.1.pep.all.fa.gz"

# Bundled DLA alleles shipped with the package
_PACKAGE_DATA_DIR = Path(__file__).parent  # dogneo/data/
_BUNDLED_DLA_ALLELES = _PACKAGE_DATA_DIR / "demo" / "dla_alleles.txt"


class ReferenceDataManager:
    """Manages canine reference data downloads and caching.

    Default cache location: ``~/.dogneo/data/``
    Override with ``DOGNEO_DATA_DIR`` env var or ``cache_dir`` parameter.

    Usage::

        mgr = ReferenceDataManager()
        proteome_path = mgr.setup()  # downloads if needed
        alleles_path = mgr.get_dla_alleles_path()
    """

    def __init__(self, cache_dir: Path | None = None) -> None:
        if cache_dir is not None:
            self.cache_dir = Path(cache_dir)
        else:
            env_dir = os.environ.get("DOGNEO_DATA_DIR")
            if env_dir:
                self.cache_dir = Path(env_dir)
            else:
                self.cache_dir = Path.home() / ".dogneo" / "data"

    # ------------------------------------------------------------------
    # Proteome
    # ------------------------------------------------------------------

    def get_proteome_path(self) -> Path | None:
        """Return cached proteome FASTA path, or None if not downloaded."""
        fa = self.cache_dir / PROTEOME_FILENAME
        return fa if fa.exists() else None

    def setup(self, force: bool = False) -> Path:
        """Download and cache CanFam3.1 proteome from Ensembl.

        Args:
            force: Re-download even if cached file exists.

        Returns:
            Path to the decompressed FASTA file.

        Raises:
            RuntimeError: If the download fails.
        """
        fa_path = self.cache_dir / PROTEOME_FILENAME
        gz_path = self.cache_dir / PROTEOME_GZ_FILENAME

        # Skip if already cached (unless force)
        if fa_path.exists() and not force:
            logger.info("Proteome already cached: %s", fa_path)
            return fa_path

        # Create cache dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Download
        logger.info("Downloading CanFam3.1 proteome from Ensembl...")
        with requests.get(PROTEOME_URL, stream=True, timeout=120) as resp:
            if not resp.ok:
                raise RuntimeError(
                    f"Failed to download proteome: {resp.status_code} {resp.text[:200]}"
                )

            total = int(resp.headers.get("content-length", 0))
            with open(gz_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)

        # Decompress
        logger.info("Decompressing %s...", gz_path.name)
        with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Clean up .gz
        gz_path.unlink(missing_ok=True)

        n_proteins = self._count_proteins(fa_path)
        logger.info("Proteome ready: %s (%d proteins)", fa_path, n_proteins)
        return fa_path

    # ------------------------------------------------------------------
    # DLA alleles
    # ------------------------------------------------------------------

    def get_dla_alleles_path(self) -> Path:
        """Return path to bundled DLA alleles file."""
        return _BUNDLED_DLA_ALLELES

    # ------------------------------------------------------------------
    # Status
    # ------------------------------------------------------------------

    def status(self) -> dict:
        """Report status of cached reference data.

        Returns:
            Dict with keys: proteome_ready, proteome_path, proteome_proteins,
            dla_alleles_ready, dla_alleles_path.
        """
        fa = self.cache_dir / PROTEOME_FILENAME
        proteome_ready = fa.exists()

        result: dict = {
            "proteome_ready": proteome_ready,
            "proteome_path": str(fa) if proteome_ready else None,
            "proteome_proteins": self._count_proteins(fa) if proteome_ready else 0,
            "dla_alleles_ready": _BUNDLED_DLA_ALLELES.exists(),
            "dla_alleles_path": str(_BUNDLED_DLA_ALLELES),
        }
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _count_proteins(fasta_path: Path) -> int:
        """Count number of protein entries (lines starting with >) in FASTA."""
        count = 0
        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
        return count
