from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import polars as pl

SORT_COLUMNS = {"rank", "risk_score", "mismatch_count", "gene_name", "transcript_id"}


@dataclass(frozen=True)
class CandidatePage:
    rows: list[dict[str, Any]]
    total: int
    page: int
    page_size: int
    pages: int
    sort: str
    direction: str


def query_candidates(
    path: Path,
    *,
    gene: str = "",
    transcript: str = "",
    regions: tuple[str, ...] = (),
    mismatches: tuple[int, ...] = (),
    min_risk: float = 0.0,
    sort: str = "rank",
    direction: str = "asc",
    page: int = 1,
    page_size: int = 25,
) -> CandidatePage:
    query = pl.scan_parquet(path)
    if gene.strip():
        needle = gene.strip().lower()
        query = query.filter(
            pl.col("gene_name").fill_null("").str.to_lowercase().str.contains(needle, literal=True)
            | pl.col("gene_id").fill_null("").str.to_lowercase().str.contains(needle, literal=True)
        )
    if transcript.strip():
        needle = transcript.strip().lower()
        query = query.filter(
            pl.col("transcript_id").fill_null("").str.to_lowercase().str.contains(needle, literal=True)
        )
    if regions:
        query = query.filter(pl.col("transcript_region").is_in(regions))
    if mismatches:
        query = query.filter(pl.col("mismatch_count").is_in(mismatches))
    query = query.filter(pl.col("risk_score") >= min_risk)

    total = int(query.select(pl.len()).collect().item())
    page_size = page_size if page_size in {25, 50, 100} else 25
    pages = max(1, math.ceil(total / page_size))
    page = min(max(page, 1), pages)
    sort = sort if sort in SORT_COLUMNS else "rank"
    direction = direction if direction in {"asc", "desc"} else "asc"
    sort_columns = [sort] if sort == "rank" else [sort, "rank"]
    descending = [direction == "desc"] + ([False] if sort != "rank" else [])
    rows = (
        query.sort(sort_columns, descending=descending, nulls_last=True)
        .slice((page - 1) * page_size, page_size)
        .collect()
        .to_dicts()
    )
    return CandidatePage(rows, total, page, page_size, pages, sort, direction)


def get_candidate(path: Path, rank: int) -> dict[str, Any]:
    rows = pl.scan_parquet(path).filter(pl.col("rank") == rank).limit(1).collect().to_dicts()
    if not rows:
        raise KeyError(f"Candidate rank {rank} was not found")
    return rows[0]


__all__ = ["CandidatePage", "get_candidate", "query_candidates"]
