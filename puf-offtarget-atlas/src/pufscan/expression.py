from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pufscan.gencode import strip_version


def tissue_specificity_tau(values: Sequence[float]) -> float | None:
    if len(values) <= 1:
        return None
    maximum = max(values)
    if maximum <= 0:
        return None
    return sum(1.0 - value / maximum for value in values) / (len(values) - 1)


def _read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith(".parquet"):
        return pd.read_parquet(path)
    separator = "\t" if name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    compression = "gzip" if name.endswith(".gz") else None
    return pd.read_csv(path, sep=separator, compression=compression)


def prepare_expression(input_path: Path, output_path: Path, id_column: str = "gene_id") -> Path:
    frame = _read_table(input_path)
    if id_column not in frame.columns:
        raise ValueError(f"Expression input lacks required ID column {id_column!r}")
    original_column = f"original_{id_column}"
    frame.insert(0, original_column, frame[id_column].astype(str))
    frame[id_column] = frame[id_column].astype(str).map(strip_version)
    duplicates = frame.loc[frame[id_column].duplicated(keep=False), id_column].unique().tolist()
    if duplicates:
        raise ValueError(f"Duplicate IDs after version stripping: {duplicates[:10]}")
    tissue_columns = [column for column in frame.columns if column not in {id_column, original_column}]
    if not tissue_columns:
        raise ValueError("Expression input contains no tissue columns")
    frame[tissue_columns] = frame[tissue_columns].apply(pd.to_numeric, errors="raise")
    values = frame[tissue_columns].to_numpy(dtype=float).ravel()
    finite = values[np.isfinite(values)]
    p99 = float(np.percentile(finite, 99)) if len(finite) else 0.0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(output_path, index=False)
    output_path.with_suffix(output_path.suffix + ".manifest.json").write_text(
        json.dumps({"id_column": id_column, "tissue_columns": tissue_columns, "p99_tpm": p99}, indent=2),
        encoding="utf-8",
    )
    return output_path


@dataclass(frozen=True)
class ExpressionSummary:
    target_tissue_tpm: float | None
    target_tissue_tpm_by_tissue: dict[str, float] | None
    all_tissue_tpm: dict[str, float] | None
    max_tissue_tpm: float | None
    median_tissue_tpm: float | None
    top_expressed_tissue: str | None
    number_of_tissues_tpm_ge_1: int | None
    number_of_tissues_tpm_ge_10: int | None
    tissue_specificity_tau: float | None
    expression_score: float | None


class ExpressionMatrix:
    def __init__(self, frame: pd.DataFrame, id_column: str, tissue_columns: list[str], p99_tpm: float):
        self.frame = frame.set_index(id_column, drop=False)
        self.id_column = id_column
        self.tissue_columns = tissue_columns
        self.p99_tpm = p99_tpm

    @classmethod
    def load(cls, path: Path) -> ExpressionMatrix:
        manifest = json.loads(path.with_suffix(path.suffix + ".manifest.json").read_text(encoding="utf-8"))
        return cls(pd.read_parquet(path), manifest["id_column"], manifest["tissue_columns"], manifest["p99_tpm"])

    def summarize(self, stable_id: str, target_tissues: Sequence[str] = ()) -> ExpressionSummary:
        stable_id = str(strip_version(stable_id))
        if stable_id not in self.frame.index:
            return ExpressionSummary(None, None, None, None, None, None, None, None, None, None)
        missing_tissues = sorted(set(target_tissues) - set(self.tissue_columns))
        if missing_tissues:
            raise ValueError(f"Unknown target tissues: {', '.join(missing_tissues)}")
        row = self.frame.loc[stable_id]
        values = [float(row[column]) for column in self.tissue_columns if not pd.isna(row[column])]
        if not values:
            return ExpressionSummary(None, None, None, None, None, None, None, None, None, None)
        by_tissue = {tissue: float(row[tissue]) for tissue in target_tissues if not pd.isna(row[tissue])}
        all_tissues = {tissue: float(row[tissue]) for tissue in self.tissue_columns if not pd.isna(row[tissue])}
        target = max(by_tissue.values()) if by_tissue else max(values)
        top_index = max(range(len(self.tissue_columns)), key=lambda index: float(row[self.tissue_columns[index]]))
        denominator = math.log1p(self.p99_tpm) if self.p99_tpm > 0 else 0.0
        score = min(math.log1p(target) / denominator, 1.0) if denominator > 0 else None
        return ExpressionSummary(
            target_tissue_tpm=target,
            target_tissue_tpm_by_tissue=by_tissue or None,
            all_tissue_tpm=all_tissues,
            max_tissue_tpm=max(values),
            median_tissue_tpm=float(np.median(values)),
            top_expressed_tissue=self.tissue_columns[top_index],
            number_of_tissues_tpm_ge_1=sum(value >= 1 for value in values),
            number_of_tissues_tpm_ge_10=sum(value >= 10 for value in values),
            tissue_specificity_tau=tissue_specificity_tau(values),
            expression_score=score,
        )
