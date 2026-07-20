from pathlib import Path

from typer.testing import CliRunner

from pufscan.cli import app

DATA = Path(__file__).parent / "data"
RUNNER = CliRunner()


def test_cli_help_lists_required_commands() -> None:
    result = RUNNER.invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in ["download-gencode", "prepare-gencode", "prepare-gtex", "scan", "report", "doctor"]:
        assert command in result.stdout


def test_cli_scan_runs_on_synthetic_data(tmp_path: Path) -> None:
    result = RUNNER.invoke(
        app,
        [
            "scan",
            "--query",
            "AACGUCUAUA",
            "--max-mismatches",
            "1",
            "--gencode-fasta",
            str(DATA / "synthetic.fa"),
            "--gencode-gtf",
            str(DATA / "synthetic.gtf"),
            "--gtex-expression",
            str(DATA / "expression.tsv"),
            "--target-tissue",
            "Liver",
            "--mode",
            "binding_only",
            "--no-structure",
            "--output",
            str(tmp_path),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "candidate sites: 6" in result.stdout


def test_editor_cli_requires_window(tmp_path: Path) -> None:
    result = RUNNER.invoke(
        app,
        [
            "scan",
            "--query",
            "AACGUCUAUA",
            "--gencode-fasta",
            str(DATA / "synthetic.fa"),
            "--gencode-gtf",
            str(DATA / "synthetic.gtf"),
            "--mode",
            "editor_fusion",
            "--editor",
            "APOBEC_C2U",
            "--output",
            str(tmp_path),
        ],
    )
    assert result.exit_code != 0

