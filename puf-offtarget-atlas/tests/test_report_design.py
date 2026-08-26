from pathlib import Path

import pandas as pd

from pufscan.config import ScanConfig, StructureConfig
from pufscan.pipeline import run_scan
from pufscan.plotting import mismatch_distribution
from pufscan.report import generate_report

DATA = Path(__file__).parent / "data"


def test_report_template_is_packaged_next_to_report_module() -> None:
    template = Path(generate_report.__code__.co_filename).parent / "templates/report.html.j2"
    assert template.is_file()


def test_offline_report_uses_atlas_design_and_embeds_dependencies(tmp_path: Path) -> None:
    result = run_scan(
        ScanConfig(
            query="AACGUCUAUA",
            fasta=DATA / "synthetic.fa",
            gtf=DATA / "synthetic.gtf",
            max_mismatches=1,
            structure=StructureConfig(enabled=False),
            output_dir=tmp_path,
            species="Testus organismus",
            genome_build="Demo1",
            annotation_provider="Test",
            annotation_release="1",
        )
    )

    html = (result.output_dir / "report.html").read_text(encoding="utf-8")

    assert 'data-report-brand="puf-atlas"' in html
    assert "Top 100 candidates" in html
    assert "Testus organismus" in html
    assert 'class="base base-a"' in html
    assert "window.Plotly" in html
    assert "<script src=" not in html


def test_plotly_figures_use_atlas_palette() -> None:
    html = mismatch_distribution(pd.DataFrame({"mismatch_count": [0, 1, 1, 2]}))

    assert "#0b7a75" in html
