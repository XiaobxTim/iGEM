from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

BRAIN_PANEL_COLUMNS = (
    "site_id",
    "gene",
    "initial_pool",
    "binding_score",
    "accessibility",
    "context_score",
    "validation_priority",
    "notes",
)

_LOCUS_COLUMNS = ("chromosome", "strand", "genomic_blocks")


def _canonical_blocks(value: object) -> str:
    if not isinstance(value, str):
        return str(value)
    try:
        return json.dumps(json.loads(value), sort_keys=True, separators=(",", ":"))
    except (TypeError, ValueError, json.JSONDecodeError):
        return value.strip()


def _site_id(row: pd.Series) -> str:
    location = "|".join(
        (
            str(row.get("chromosome", "")),
            str(row.get("strand", "")),
            _canonical_blocks(row.get("genomic_blocks", "")),
        )
    )
    digest = hashlib.sha256(location.encode("utf-8")).hexdigest()[:12]
    return f"puf-{digest}"


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    return bool(pd.isna(value))


def _score(value: object, fallback: float) -> float:
    if _is_missing(value):
        return fallback
    return max(0.0, min(1.0, float(str(value))))


def _metadata(
    *,
    source_count: int,
    unique_count: int,
    exported_count: int,
    limit: int,
    expression_fallback_ids: list[str],
    accessibility_fallback_ids: list[str],
) -> dict[str, Any]:
    warnings = []
    if expression_fallback_ids:
        warnings.append(
            "Missing expression scores were conservatively replaced with initial_pool=1.0."
        )
    if accessibility_fallback_ids:
        warnings.append(
            "Missing accessibility scores were conservatively replaced with accessibility=1.0."
        )
    return {
        "format": "puf-to-brain-candidate-panel",
        "format_version": "1.0",
        "generated_at": datetime.now(UTC).isoformat(),
        "source_candidate_count": source_count,
        "unique_locus_count": unique_count,
        "exported_locus_count": exported_count,
        "limit": limit,
        "mapping": {
            "site_id": "stable hash of chromosome, strand, and genomic_blocks",
            "gene": "gene_name, falling back to gene_id",
            "initial_pool": "expression_score",
            "binding_score": "sequence_score",
            "accessibility": "accessibility_score",
            "context_score": "consequence_score",
            "validation_priority": "risk_priority",
        },
        "conservative_fallbacks": {
            "expression_score": {
                "replacement": 1.0,
                "count": len(expression_fallback_ids),
                "site_ids": expression_fallback_ids,
            },
            "accessibility_score": {
                "replacement": 1.0,
                "count": len(accessibility_fallback_ids),
                "site_ids": accessibility_fallback_ids,
            },
        },
        "warnings": warnings,
        "disclaimer": (
            "This panel is a prioritization bridge between two literature-informed models. "
            "It is not a calibrated probability, clinical prediction, or substitute for validation."
        ),
    }


def export_brain_candidate_panel(
    frame: pd.DataFrame,
    run_dir: Path,
    limit: int = 100,
) -> tuple[Path, Path]:
    """Write the highest-risk unique genomic loci in Brain model CSV format."""
    if limit < 1:
        raise ValueError("limit must be at least 1")

    csv_path = run_dir / "brain_candidate_panel.csv"
    metadata_path = run_dir / "brain_candidate_panel.metadata.json"
    panel_rows: list[dict[str, object]] = []
    expression_fallback_ids: list[str] = []
    accessibility_fallback_ids: list[str] = []

    if frame.empty:
        unique_count = 0
        selected = frame
    else:
        missing_columns = set(_LOCUS_COLUMNS) - set(frame.columns)
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"candidate frame is missing locus columns: {missing}")
        selected = (
            frame.sort_values(["risk_score", "sequence_score"], ascending=False)
            .drop_duplicates(list(_LOCUS_COLUMNS), keep="first")
        )
        unique_count = len(selected)
        selected = selected.head(limit)

    for _, row in selected.iterrows():
        site_id = _site_id(row)
        expression_missing = _is_missing(row.get("expression_score"))
        accessibility_missing = _is_missing(row.get("accessibility_score"))
        if expression_missing:
            expression_fallback_ids.append(site_id)
        if accessibility_missing:
            accessibility_fallback_ids.append(site_id)

        gene_name = row.get("gene_name")
        gene = row.get("gene_id", "") if _is_missing(gene_name) else gene_name
        notes = [
            f"PUF transcript {row.get('transcript_id', '')}",
            f"risk score {float(row.get('risk_score', 0.0)):.2f}/100",
        ]
        if expression_missing:
            notes.append("conservative expression fallback")
        if accessibility_missing:
            notes.append("conservative accessibility fallback")
        panel_rows.append(
            {
                "site_id": site_id,
                "gene": str(gene),
                "initial_pool": _score(row.get("expression_score"), 1.0),
                "binding_score": _score(row.get("sequence_score"), 0.0),
                "accessibility": _score(row.get("accessibility_score"), 1.0),
                "context_score": _score(row.get("consequence_score"), 0.0),
                "validation_priority": str(row.get("risk_priority", "Unclassified")),
                "notes": "; ".join(notes),
            }
        )

    panel = pd.DataFrame(panel_rows, columns=BRAIN_PANEL_COLUMNS)
    panel.to_csv(csv_path, index=False)
    metadata = _metadata(
        source_count=len(frame),
        unique_count=unique_count,
        exported_count=len(panel),
        limit=limit,
        expression_fallback_ids=expression_fallback_ids,
        accessibility_fallback_ids=accessibility_fallback_ids,
    )
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return csv_path, metadata_path
