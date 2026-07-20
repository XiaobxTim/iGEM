from __future__ import annotations

import csv
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from pufscan.sequence import hamming_distance


@dataclass(frozen=True)
class SequenceScore:
    exact_match: bool
    mismatch_count: int
    mismatch_positions_1based: tuple[int, ...]
    weighted_mismatch_penalty: float
    sequence_identity: float
    sequence_score: float


class SequenceScoringModel(Protocol):
    def score(self, target: str, observed: str) -> SequenceScore: ...


@dataclass(frozen=True)
class UniformMismatchScorer:
    weights: tuple[float, ...] | None = None

    def score(self, target: str, observed: str) -> SequenceScore:
        mismatch_count, positions = hamming_distance(target, observed)
        weights = self.weights or (1.0,) * len(target)
        if len(weights) != len(target):
            raise ValueError("Position weight count must equal motif length")
        if any(weight <= 0 for weight in weights):
            raise ValueError("Position weights must be positive")
        penalty = sum(weights[position - 1] for position in positions)
        return SequenceScore(
            exact_match=mismatch_count == 0,
            mismatch_count=mismatch_count,
            mismatch_positions_1based=positions,
            weighted_mismatch_penalty=penalty,
            sequence_identity=1.0 - mismatch_count / len(target),
            sequence_score=max(0.0, min(1.0, 1.0 - penalty / sum(weights))),
        )


@dataclass(frozen=True)
class SubstitutionMatrixScorer:
    compatibility: Mapping[tuple[int, str, str], float]
    weights: tuple[float, ...] | None = None

    def score(self, target: str, observed: str) -> SequenceScore:
        mismatch_count, positions = hamming_distance(target, observed)
        weights = self.weights or (1.0,) * len(target)
        if len(weights) != len(target):
            raise ValueError("Position weight count must equal motif length")
        values: list[float] = []
        for position, (expected, actual) in enumerate(zip(target, observed, strict=True), 1):
            key = (position, expected, actual)
            if key not in self.compatibility:
                raise ValueError(f"Substitution matrix lacks entry {key}")
            values.append(float(self.compatibility[key]))
        score = sum(weight * value for weight, value in zip(weights, values, strict=True)) / sum(weights)
        penalty = sum(weights[position - 1] for position in positions)
        return SequenceScore(
            exact_match=mismatch_count == 0,
            mismatch_count=mismatch_count,
            mismatch_positions_1based=positions,
            weighted_mismatch_penalty=penalty,
            sequence_identity=1.0 - mismatch_count / len(target),
            sequence_score=max(0.0, min(1.0, score)),
        )


def load_position_weights(path: Path, motif_length: int) -> tuple[float, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    weights = {int(row["position"]): float(row["weight"]) for row in rows}
    if set(weights) != set(range(1, motif_length + 1)):
        raise ValueError("Position weight file must contain every motif position exactly once")
    return tuple(weights[position] for position in range(1, motif_length + 1))


def load_substitution_matrix(
    path: Path, motif_length: int
) -> dict[tuple[int, str, str], float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    matrix = {
        (int(row["position"]), row["target_base"].upper(), row["observed_base"].upper()): float(
            row["compatibility_score"]
        )
        for row in rows
    }
    expected = {
        (position, target, observed)
        for position in range(1, motif_length + 1)
        for target in "ACGT"
        for observed in "ACGT"
    }
    missing = expected - set(matrix)
    extra = set(matrix) - expected
    if missing or extra:
        raise ValueError(
            f"Substitution matrix must contain all 16 base pairs at every position; "
            f"missing={len(missing)}, extra={len(extra)}"
        )
    if any(not 0 <= score <= 1 for score in matrix.values()):
        raise ValueError("Substitution compatibility scores must be within [0, 1]")
    return matrix


@dataclass(frozen=True)
class RiskScoreResult:
    risk_score: float
    risk_priority: str
    missing_features: tuple[str, ...]
    normalized_weights: dict[str, float]


def risk_priority(score: float) -> str:
    if score < 25:
        return "Low priority"
    if score < 50:
        return "Moderate priority"
    if score < 75:
        return "High priority"
    return "Very high priority"


def calculate_risk_score(
    components: Mapping[str, float | None], weights: Mapping[str, float]
) -> RiskScoreResult:
    missing = tuple(name for name in weights if components.get(name) is None)
    present = {name: weight for name, weight in weights.items() if components.get(name) is not None and weight > 0}
    total_weight = sum(present.values())
    if total_weight <= 0:
        raise ValueError("At least one weighted risk component must be available")
    normalized = {name: weight / total_weight for name, weight in present.items()}
    log_fraction = 0.0
    for name, normalized_weight in normalized.items():
        component = components[name]
        if component is None:
            raise AssertionError("Normalized risk component cannot be missing")
        log_fraction += normalized_weight * math.log(max(0.001, min(1.0, component)))
    score = 100.0 * math.exp(log_fraction)
    return RiskScoreResult(score, risk_priority(score), missing, normalized)
