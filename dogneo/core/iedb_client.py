"""IEDB API client for MHC/DLA binding prediction.

Free, no-install binding prediction via the IEDB web API.
Supports DLA alleles through the NetMHCpan EL backend.
Includes disk caching and rate limiting for responsible usage.
"""
from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path

import requests

from dogneo.core.binding import BindingPrediction

logger = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path.home() / ".dogneo" / "cache" / "iedb"
IEDB_MHCI_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
REQUEST_INTERVAL = 1.0  # seconds between API calls


class IEDBClient:
    """IEDB MHC-I binding prediction client with caching.

    Args:
        cache_dir: Directory for caching prediction results.
        method: IEDB prediction method (default: recommended).
        timeout: HTTP request timeout in seconds.
    """

    def __init__(
        self,
        cache_dir: Path | None = None,
        method: str = "recommended",
        timeout: int = 120,
    ):
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        self.method = method
        self.timeout = timeout
        self._last_request_time: float = 0.0

    def _cache_key(self, peptides: list[str], alleles: list[str]) -> str:
        """Generate a deterministic cache key from inputs."""
        raw = f"{sorted(peptides)}|{sorted(alleles)}|{self.method}"
        return hashlib.sha256(raw.encode()).hexdigest()[:16]

    def _load_cache(self, key: str) -> list[BindingPrediction] | None:
        """Load cached predictions if available."""
        cache_file = self.cache_dir / f"{key}.json"
        if not cache_file.exists():
            return None

        try:
            data = json.loads(cache_file.read_text())
            return [
                BindingPrediction(
                    peptide_sequence=d["peptide_sequence"],
                    allele=d["allele"],
                    affinity_nm=d["affinity_nm"],
                    percentile_rank=d["percentile_rank"],
                    tool="iedb",
                    mhc_class=d.get("mhc_class", 1),
                )
                for d in data
            ]
        except (json.JSONDecodeError, KeyError) as e:
            logger.debug("Cache read failed for %s: %s", key, e)
            return None

    def _save_cache(self, key: str, predictions: list[BindingPrediction]) -> None:
        """Save predictions to disk cache."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = self.cache_dir / f"{key}.json"
        data = [
            {
                "peptide_sequence": p.peptide_sequence,
                "allele": p.allele,
                "affinity_nm": p.affinity_nm,
                "percentile_rank": p.percentile_rank,
                "mhc_class": p.mhc_class,
            }
            for p in predictions
        ]
        cache_file.write_text(json.dumps(data))

    def _rate_limit(self) -> None:
        """Enforce minimum interval between API requests."""
        elapsed = time.monotonic() - self._last_request_time
        if elapsed < REQUEST_INTERVAL:
            time.sleep(REQUEST_INTERVAL - elapsed)

    def predict_batch(
        self,
        peptides: list[str],
        alleles: list[str],
    ) -> list[BindingPrediction]:
        """Predict MHC-I binding for peptides against alleles via IEDB API.

        Args:
            peptides: Peptide sequences to evaluate.
            alleles: MHC/DLA allele names.

        Returns:
            List of BindingPrediction results.

        Raises:
            requests.ConnectionError: If IEDB API is unreachable.
        """
        if not peptides or not alleles:
            return []

        # Check cache
        cache_key = self._cache_key(peptides, alleles)
        cached = self._load_cache(cache_key)
        if cached is not None:
            logger.info("IEDB cache hit (%d predictions)", len(cached))
            return cached

        # Rate limit
        self._rate_limit()

        # Build request
        sequences = "\n".join(peptides)
        allele_str = ",".join(alleles)
        length_set = sorted(set(len(p) for p in peptides))
        length_str = ",".join(str(n) for n in length_set)

        response = requests.post(
            IEDB_MHCI_URL,
            data={
                "method": self.method,
                "sequence_text": sequences,
                "allele": allele_str,
                "length": length_str,
            },
            timeout=self.timeout,
        )
        self._last_request_time = time.monotonic()

        predictions: list[BindingPrediction] = []

        if response.ok:
            lines = response.text.strip().split("\n")
            if len(lines) > 1:
                header = lines[0].split("\t")
                for line in lines[1:]:
                    fields = line.split("\t")
                    if len(fields) < len(header):
                        continue
                    row = dict(zip(header, fields))  # noqa: B905
                    try:
                        predictions.append(
                            BindingPrediction(
                                peptide_sequence=row.get("peptide", ""),
                                allele=row.get("allele", ""),
                                affinity_nm=float(row.get("ic50", "99999")),
                                percentile_rank=float(row.get("percentile_rank", "100")),
                                tool="iedb",
                                mhc_class=1,
                            )
                        )
                    except (ValueError, KeyError):
                        continue
        else:
            logger.error("IEDB API error: %d %s", response.status_code, response.text[:200])

        # Cache results
        if predictions:
            self._save_cache(cache_key, predictions)

        logger.info("IEDB returned %d predictions", len(predictions))
        return predictions
