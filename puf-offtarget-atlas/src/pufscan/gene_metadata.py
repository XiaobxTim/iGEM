from __future__ import annotations

import json
import time
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
import requests

from pufscan.gencode import strip_version


@dataclass(frozen=True)
class GeneMetadata:
    gene_id: str
    gene_name: str | None = None
    description: str = "Not available"
    go_terms: str = "Not available"
    pathways: str = "Not available"
    disease_annotations: str = "Not available"
    essentiality: str = "Not available"


class GeneMetadataProvider(ABC):
    @abstractmethod
    def get_gene_metadata(self, gene_id: str) -> GeneMetadata: ...


class LocalGeneMetadataProvider(GeneMetadataProvider):
    def __init__(self, path: Path):
        separator = "\t" if path.name.endswith((".tsv", ".tsv.gz")) else ","
        frame = pd.read_csv(path, sep=separator)
        if "gene_id" not in frame.columns:
            raise ValueError("Gene metadata file lacks gene_id")
        frame["stable_gene_id"] = frame["gene_id"].astype(str).map(strip_version)
        if frame["stable_gene_id"].duplicated().any():
            raise ValueError("Gene metadata contains duplicate stable gene IDs")
        self.frame = frame.set_index("stable_gene_id")

    def get_gene_metadata(self, gene_id: str) -> GeneMetadata:
        stable = str(strip_version(gene_id))
        if stable not in self.frame.index:
            return GeneMetadata(gene_id=gene_id)
        row = self.frame.loc[stable]
        return GeneMetadata(
            gene_id=gene_id,
            gene_name=_available(row.get("gene_name"), None),
            description=_available_str(row.get("description")),
            go_terms=_available_str(row.get("go_terms")),
            pathways=_available_str(row.get("pathways")),
            disease_annotations=_available_str(row.get("disease_annotations")),
            essentiality=_available_str(row.get("essentiality")),
        )


def _available(value: object, fallback: str | None = "Not available") -> str | None:
    return fallback if value is None or pd.isna(value) or str(value).strip() == "" else str(value)


def _available_str(value: object) -> str:
    return str(_available(value, "Not available"))


class EnsemblRestGeneMetadataProvider(GeneMetadataProvider):
    def __init__(self, cache_dir: Path, timeout: float = 10.0, retries: int = 3):
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.retries = retries

    def get_gene_metadata(self, gene_id: str) -> GeneMetadata:
        stable = str(strip_version(gene_id))
        cache_path = self.cache_dir / f"{stable}.json"
        if cache_path.exists():
            return GeneMetadata(**json.loads(cache_path.read_text(encoding="utf-8")))
        for attempt in range(self.retries):
            try:
                response = requests.get(
                    f"https://rest.ensembl.org/lookup/id/{stable}",
                    headers={"Content-Type": "application/json"},
                    timeout=self.timeout,
                )
                response.raise_for_status()
                payload = response.json()
                result = GeneMetadata(
                    gene_id=gene_id,
                    gene_name=payload.get("display_name"),
                    description=payload.get("description") or "Not available",
                )
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
                return result
            except (requests.RequestException, ValueError):
                if attempt + 1 < self.retries:
                    time.sleep(0.5 * 2**attempt)
        return GeneMetadata(gene_id=gene_id)
