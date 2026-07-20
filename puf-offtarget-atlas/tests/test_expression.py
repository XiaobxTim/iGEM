from pathlib import Path

import pytest

from pufscan.expression import ExpressionMatrix, prepare_expression, tissue_specificity_tau

DATA = Path(__file__).parent / "data"


def test_tau_handles_regular_and_all_zero_expression() -> None:
    assert tissue_specificity_tau([100.0, 10.0, 1.0]) == pytest.approx(0.945)
    assert tissue_specificity_tau([0.0, 0.0, 0.0]) is None
    assert tissue_specificity_tau([5.0]) is None


def test_prepare_expression_strips_versions_and_preserves_na(tmp_path: Path) -> None:
    output = tmp_path / "expression.parquet"
    prepare_expression(DATA / "expression.tsv", output)
    matrix = ExpressionMatrix.load(output)
    summary = matrix.summarize("ENSG000001", ("Liver",))
    assert summary.target_tissue_tpm == pytest.approx(100.0)
    assert summary.top_expressed_tissue == "Liver"
    assert summary.number_of_tissues_tpm_ge_10 == 2
    assert summary.expression_score is not None
    assert matrix.summarize("ENSG_MISSING", ("Liver",)).max_tissue_tpm is None

