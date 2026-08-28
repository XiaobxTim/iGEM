from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRAIN_ROOT = ROOT / "brain_delivery_model" / "XiaobxTim-iGEM-REWIRE-iterated-v4"
PUF_ROOT = ROOT / "puf-offtarget-atlas"
sys.path.insert(0, str(PUF_ROOT / "src"))
sys.path.insert(0, str(BRAIN_ROOT))

from brain_app.service import parse_candidate_panel, run_model  # noqa: E402
from pufscan.config import ScanConfig, StructureConfig  # noqa: E402
from pufscan.pipeline import run_scan  # noqa: E402


def test_synthetic_puf_panel_runs_in_brain_model(tmp_path: Path) -> None:
    data = PUF_ROOT / "tests" / "data"
    result = run_scan(
        ScanConfig(
            query="AACGUCUAUA",
            fasta=data / "synthetic.fa",
            gtf=data / "synthetic.gtf",
            expression=data / "expression.tsv",
            target_tissues=("Liver",),
            max_mismatches=1,
            structure=StructureConfig(enabled=False),
            output_dir=tmp_path,
        )
    )
    panel_path = result.output_dir / "brain_candidate_panel.csv"

    rows = parse_candidate_panel(panel_path.read_bytes())
    model_result = run_model(
        mode="single",
        design_id="10R-design-4",
        route="footpad",
        dose=1.0,
        duration=24.0,
        panel_rows=rows,
    )

    assert rows
    assert model_result["panel_summary"]["n_sites"] == float(len(rows))
    assert model_result["panel_summary"]["provided"] is True
    assert 0.0 <= model_result["metrics"]["apoe3_like_fraction_final"] <= 1.0
    assert model_result["metrics"]["off_target_burden_final"] >= 0.0
