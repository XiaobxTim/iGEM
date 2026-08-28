from __future__ import annotations

import csv

import numpy as np
import pytest

from spatial_model.pbpk_boundary import (
    BoundaryCurve,
    build_boundary_curves,
    simulate_module12,
    write_boundary_csv,
)
from utils.config_loader import load_base_config


def test_module12_curve_is_nonnegative_and_dose_ordered():
    config = load_base_config(".")

    low = simulate_module12(config, dose=0.3, t_end_h=12.0, dt_h=0.5)
    high = simulate_module12(config, dose=3.0, t_end_h=12.0, dt_h=0.5)

    assert np.all(low.time_min >= 0.0)
    assert np.all(np.diff(low.time_min) > 0.0)
    assert np.all(low.raw_amount >= 0.0)
    assert np.all(high.raw_amount >= low.raw_amount)
    assert high.raw_amount.max() > low.raw_amount.max()


def test_three_curves_share_medium_dose_normalization():
    config = load_base_config(".")

    curves = build_boundary_curves(
        config,
        doses=(0.3, 1.0, 3.0),
        t_end_h=12.0,
        dt_h=0.5,
    )

    assert set(curves) == {0.3, 1.0, 3.0}
    assert curves[1.0].normalized_amount.max() == pytest.approx(1.0)
    for curve in curves.values():
        assert np.allclose(
            curve.normalized_amount,
            curve.raw_amount / curves[1.0].raw_amount.max(),
        )


def test_boundary_csv_has_stable_schema(tmp_path):
    curve = BoundaryCurve(
        dose=1.0,
        time_min=np.array([0.0, 30.0, 60.0]),
        raw_amount=np.array([0.0, 0.2, 0.1]),
        normalized_amount=np.array([0.0, 1.0, 0.5]),
    )
    output = tmp_path / "brain_blood_input.csv"

    write_boundary_csv(curve, output)

    assert b"\r" not in output.read_bytes()
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "time_min",
        "A_brain_blood_raw",
        "A_brain_blood_normalized",
    ]
    assert rows[-1] == {
        "time_min": "60",
        "A_brain_blood_raw": "0.1",
        "A_brain_blood_normalized": "0.5",
    }


@pytest.mark.parametrize(
    ("dose", "t_end_h", "dt_h"),
    [(-1.0, 12.0, 0.5), (1.0, 0.0, 0.5), (1.0, 12.0, 0.0)],
)
def test_module12_rejects_invalid_run_arguments(dose, t_end_h, dt_h):
    config = load_base_config(".")

    with pytest.raises(ValueError):
        simulate_module12(config, dose=dose, t_end_h=t_end_h, dt_h=dt_h)
