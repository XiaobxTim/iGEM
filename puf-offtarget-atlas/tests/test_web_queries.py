from pathlib import Path

import pandas as pd

from pufscan.web.queries import get_candidate, query_candidates


def make_candidates(path: Path) -> Path:
    pd.DataFrame(
        [
            {"rank": 1, "gene_name": "APOE", "gene_id": "ENSG1", "transcript_id": "ENST1", "mismatch_count": 0, "transcript_region": "CDS", "risk_score": 92.0},
            {"rank": 2, "gene_name": "APOB", "gene_id": "ENSG2", "transcript_id": "ENST2", "mismatch_count": 1, "transcript_region": "3UTR", "risk_score": 66.0},
            {"rank": 3, "gene_name": "LDLR", "gene_id": "ENSG3", "transcript_id": "ENST3", "mismatch_count": 2, "transcript_region": "CDS", "risk_score": 41.0},
        ]
    ).to_parquet(path, index=False)
    return path


def test_candidate_query_filters_and_paginates(tmp_path: Path) -> None:
    path = make_candidates(tmp_path / "hits.parquet")

    page = query_candidates(
        path,
        gene="apo",
        regions=("CDS",),
        min_risk=80,
        page=1,
        page_size=25,
    )

    assert page.total == 1
    assert page.rows[0]["gene_name"] == "APOE"
    assert page.pages == 1


def test_candidate_query_uses_stable_server_side_sort(tmp_path: Path) -> None:
    path = make_candidates(tmp_path / "hits.parquet")

    page = query_candidates(path, sort="risk_score", direction="asc", page=1, page_size=25)

    assert [row["rank"] for row in page.rows] == [3, 2, 1]
    assert page.total == 3
    assert page.pages == 1


def test_candidate_detail_is_selected_by_stable_rank(tmp_path: Path) -> None:
    path = make_candidates(tmp_path / "hits.parquet")

    candidate = get_candidate(path, 2)

    assert candidate["transcript_id"] == "ENST2"
