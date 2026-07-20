import math

import pytest

from pufscan.scoring import (
    SubstitutionMatrixScorer,
    UniformMismatchScorer,
    calculate_risk_score,
    load_substitution_matrix,
    risk_priority,
)


def test_uniform_sequence_score_supports_position_weights() -> None:
    scorer = UniformMismatchScorer(weights=(1.0, 2.0, 1.0, 1.0))
    result = scorer.score("ACGT", "ATGT")
    assert result.exact_match is False
    assert result.mismatch_count == 1
    assert result.weighted_mismatch_penalty == pytest.approx(2.0)
    assert result.sequence_identity == pytest.approx(0.75)
    assert result.sequence_score == pytest.approx(0.6)


def test_risk_score_drops_missing_weights_and_renormalizes() -> None:
    result = calculate_risk_score(
        components={"sequence": 0.8, "accessibility": None, "expression": 0.5, "consequence": 1.0},
        weights={"sequence": 0.5, "accessibility": 0.2, "expression": 0.2, "consequence": 0.1},
    )
    expected = 100 * math.prod((0.8 ** 0.625, 0.5 ** 0.25, 1.0 ** 0.125))
    assert result.risk_score == pytest.approx(expected)
    assert result.missing_features == ("accessibility",)


@pytest.mark.parametrize(
    ("score", "label"),
    [(0, "Low priority"), (25, "Moderate priority"), (50, "High priority"), (75, "Very high priority"), (100, "Very high priority")],
)
def test_risk_priority_boundaries(score: float, label: str) -> None:
    assert risk_priority(score) == label


def test_complete_substitution_matrix_can_be_loaded(tmp_path) -> None:
    matrix = tmp_path / "matrix.csv"
    rows = ["position,target_base,observed_base,compatibility_score"]
    for position in range(1, 9):
        for target in "ACGT":
            for observed in "ACGT":
                rows.append(f"{position},{target},{observed},{1.0 if target == observed else 0.2}")
    matrix.write_text("\n".join(rows), encoding="utf-8")
    values = load_substitution_matrix(matrix, 8)
    result = SubstitutionMatrixScorer(values).score("AACGTCTA", "ATCGTCTA")
    assert result.sequence_score == pytest.approx(0.9)
