from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.integrate import solve_ivp

from models.full_model.state_vector import IDX, STATE_ORDER
from models.pbpk.lymphatic_absorption import lymphatic_absorption_rhs
from models.pbpk.organ_distribution import organ_distribution_rhs


@dataclass(frozen=True)
class BoundaryCurve:
    """Prescribed brain-blood AAV input for one normalized dose."""

    dose: float
    time_min: np.ndarray
    raw_amount: np.ndarray
    normalized_amount: np.ndarray


def _validate_run_arguments(dose: float, t_end_h: float, dt_h: float) -> None:
    if not np.isfinite(dose) or dose < 0.0:
        raise ValueError("dose must be finite and non-negative")
    if not np.isfinite(t_end_h) or t_end_h <= 0.0:
        raise ValueError("t_end_h must be finite and positive")
    if not np.isfinite(dt_h) or dt_h <= 0.0 or dt_h > t_end_h:
        raise ValueError("dt_h must be finite, positive, and no larger than t_end_h")


def simulate_module12(config: dict, dose: float, t_end_h: float, dt_h: float) -> BoundaryCurve:
    """Integrate absorption and distribution only, leaving BBB transport spatial."""

    _validate_run_arguments(dose, t_end_h, dt_h)
    route = str(config["route"])
    if route not in config["absorption"]:
        raise ValueError(f"unsupported administration route: {route}")

    y0 = np.zeros(len(STATE_ORDER), dtype=float)
    if route == "iv":
        y0[IDX["A_blood"]] = dose
    else:
        y0[IDX["A_dep"]] = dose

    absorption = config["absorption"][route]
    distribution = config["distribution"]

    def rhs(time_h: float, state: np.ndarray) -> np.ndarray:
        return (
            lymphatic_absorption_rhs(time_h, state, absorption, IDX)
            + organ_distribution_rhs(time_h, state, distribution, IDX)
        )

    sample_count = int(np.floor(t_end_h / dt_h))
    t_eval = np.arange(sample_count + 1, dtype=float) * dt_h
    if t_eval[-1] < t_end_h:
        t_eval = np.append(t_eval, t_end_h)

    solution = solve_ivp(
        rhs,
        (0.0, t_end_h),
        y0,
        t_eval=t_eval,
        method=str(config.get("simulation", {}).get("solver", "LSODA")),
        rtol=1e-8,
        atol=1e-11,
    )
    if not solution.success:
        raise RuntimeError(f"Module 1-2 boundary simulation failed: {solution.message}")

    raw = np.maximum(solution.y[IDX["A_brain_blood"]], 0.0)
    if not np.all(np.isfinite(raw)):
        raise RuntimeError("Module 1-2 boundary simulation produced non-finite values")
    return BoundaryCurve(
        dose=float(dose),
        time_min=solution.t * 60.0,
        raw_amount=raw,
        normalized_amount=np.zeros_like(raw),
    )


def build_boundary_curves(
    config: dict,
    doses: Iterable[float],
    t_end_h: float,
    dt_h: float,
    reference_dose: float = 1.0,
) -> dict[float, BoundaryCurve]:
    """Simulate dose curves and normalize all of them to the medium-dose peak."""

    dose_values = tuple(float(value) for value in doses)
    if len(set(dose_values)) != len(dose_values):
        raise ValueError("doses must be unique")
    if reference_dose not in dose_values:
        raise ValueError(f"reference dose {reference_dose:g} must be included")

    curves = {
        dose: simulate_module12(config, dose=dose, t_end_h=t_end_h, dt_h=dt_h)
        for dose in dose_values
    }
    reference_peak = float(curves[reference_dose].raw_amount.max())
    if not np.isfinite(reference_peak) or reference_peak <= 0.0:
        raise RuntimeError("reference dose has no positive brain-blood AAV peak")

    return {
        dose: replace(curve, normalized_amount=curve.raw_amount / reference_peak)
        for dose, curve in curves.items()
    }


def _format_number(value: float) -> str:
    return format(float(value), ".12g")


def write_boundary_csv(curve: BoundaryCurve, path: str | Path) -> Path:
    """Write one validated boundary curve using the stable public CSV schema."""

    arrays = (curve.time_min, curve.raw_amount, curve.normalized_amount)
    if len({len(array) for array in arrays}) != 1 or not len(curve.time_min):
        raise ValueError("boundary arrays must be non-empty and have equal length")
    if not np.all(np.isfinite(np.concatenate(arrays))):
        raise ValueError("boundary curve contains non-finite values")
    if curve.time_min[0] < 0.0 or np.any(np.diff(curve.time_min) <= 0.0):
        raise ValueError("boundary time must be non-negative and strictly increasing")
    if np.any(curve.raw_amount < 0.0) or np.any(curve.normalized_amount < 0.0):
        raise ValueError("boundary AAV values must be non-negative")

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("time_min", "A_brain_blood_raw", "A_brain_blood_normalized")
        )
        for values in zip(*arrays, strict=True):
            writer.writerow(tuple(_format_number(value) for value in values))
    return output

