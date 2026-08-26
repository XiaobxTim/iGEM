import json
from pathlib import Path

import pandas as pd

from pufscan.brain_panel import BRAIN_PANEL_COLUMNS, export_brain_candidate_panel


def _candidate(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "chromosome": "chr19",
        "strand": "+",
        "genomic_blocks": '[{"chromosome":"chr19","start_1based":100,"end_1based":109}]',
        "gene_name": "APOE",
        "gene_id": "ENSG00000130203.10",
        "transcript_id": "ENST00000252486.9",
        "risk_score": 82.5,
        "risk_priority": "High priority",
        "sequence_score": 0.91,
        "expression_score": 0.64,
        "accessibility_score": 0.72,
        "consequence_score": 0.85,
    }
    row.update(overrides)
    return row


def test_export_keeps_highest_risk_row_per_locus_and_sorts_by_risk(tmp_path: Path) -> None:
    duplicate_lower_risk = _candidate(
        transcript_id="ENST_LOW",
        risk_score=30.0,
        sequence_score=0.5,
    )
    second_locus = _candidate(
        chromosome="chr2",
        genomic_blocks='[{"chromosome":"chr2","start_1based":20,"end_1based":29}]',
        gene_name="GENE2",
        gene_id="ENSG2",
        risk_score=70.0,
        risk_priority="Medium priority",
    )

    csv_path, metadata_path = export_brain_candidate_panel(
        pd.DataFrame([duplicate_lower_risk, second_locus, _candidate()]),
        tmp_path,
    )

    panel = pd.read_csv(csv_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert list(panel.columns) == list(BRAIN_PANEL_COLUMNS)
    assert list(panel["gene"]) == ["APOE", "GENE2"]
    assert panel.loc[0, "binding_score"] == 0.91
    assert panel["site_id"].is_unique
    assert panel.loc[0, "site_id"].startswith("puf-")
    assert metadata["source_candidate_count"] == 3
    assert metadata["unique_locus_count"] == 2
    assert metadata["exported_locus_count"] == 2


def test_export_uses_conservative_fallbacks_and_records_warnings(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        [
            _candidate(
                gene_name=None,
                expression_score=None,
                accessibility_score=float("nan"),
            )
        ]
    )

    csv_path, metadata_path = export_brain_candidate_panel(frame, tmp_path)

    panel = pd.read_csv(csv_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert panel.loc[0, "gene"] == "ENSG00000130203.10"
    assert panel.loc[0, "initial_pool"] == 1.0
    assert panel.loc[0, "accessibility"] == 1.0
    assert metadata["conservative_fallbacks"]["expression_score"]["count"] == 1
    assert metadata["conservative_fallbacks"]["accessibility_score"]["count"] == 1
    assert len(metadata["warnings"]) == 2


def test_export_limits_panel_and_writes_header_only_for_zero_hits(tmp_path: Path) -> None:
    rows = [
        _candidate(
            chromosome=f"chr{index}",
            genomic_blocks=(
                f'[{json.dumps({"chromosome": f"chr{index}", "start_1based": index, "end_1based": index + 9})}]'
            ),
            risk_score=float(index),
        )
        for index in range(1, 105)
    ]
    limited_dir = tmp_path / "limited"
    empty_dir = tmp_path / "empty"
    limited_dir.mkdir()
    empty_dir.mkdir()

    limited_csv, _ = export_brain_candidate_panel(pd.DataFrame(rows), limited_dir)
    empty_csv, empty_metadata = export_brain_candidate_panel(pd.DataFrame(), empty_dir)

    assert len(pd.read_csv(limited_csv)) == 100
    empty_panel = pd.read_csv(empty_csv)
    assert empty_panel.empty
    assert list(empty_panel.columns) == list(BRAIN_PANEL_COLUMNS)
    assert json.loads(empty_metadata.read_text(encoding="utf-8"))["exported_locus_count"] == 0
