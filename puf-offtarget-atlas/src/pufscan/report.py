from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from plotly.offline import get_plotlyjs

from pufscan.plotting import (
    accessibility_plot,
    mismatch_distribution,
    risk_distribution,
    tissue_heatmap,
    transcript_schematic,
)

LIMITATIONS = (
    "Sequence similarity does not prove PUF binding.",
    "Predicted accessibility does not represent in-cell RNA structure.",
    "Bulk tissue expression does not guarantee co-localization with the PUF editor.",
    "Predicted editing consequences require experimental validation.",
    "The integrated risk score is a prioritization score, not a calibrated probability.",
)


def generate_report(run_dir: Path) -> Path:
    frame = pd.read_parquet(run_dir / "all_transcript_hits.parquet")
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    template_dir = Path(__file__).resolve().parents[2] / "templates"
    environment = Environment(loader=FileSystemLoader(template_dir), autoescape=select_autoescape())
    template = environment.get_template("report.html.j2")
    html = template.render(
        summary=summary,
        plotly_js=get_plotlyjs(),
        mismatch_plot=mismatch_distribution(frame) if not frame.empty else "<p>Not available</p>",
        risk_plot=risk_distribution(frame) if not frame.empty else "<p>Not available</p>",
        schematic=transcript_schematic(frame.iloc[0]) if not frame.empty else "<p>Not available</p>",
        tissue_heatmap=tissue_heatmap(frame) if not frame.empty else "<p>Not available</p>",
        accessibility_plot=accessibility_plot(frame.iloc[0]) if not frame.empty else "<p>Not available</p>",
        top_hits=frame.head(50).to_dict(orient="records"),
        limitations=LIMITATIONS,
    )
    output = run_dir / "report.html"
    output.write_text(html, encoding="utf-8")
    return output
