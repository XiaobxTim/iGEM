from __future__ import annotations

from typing import Dict, Any, List
import copy
import numpy as np

from models.full_model.simulator import run_simulation
from optimization.metrics import extract_summary_metrics
from optimization.therapeutic_window import evaluate_therapeutic_window


def build_dose_grid(config: Dict) -> np.ndarray:
    scan_cfg = config["optimization"]["dose_scan"]
    return np.linspace(
        float(scan_cfg["dose_min"]),
        float(scan_cfg["dose_max"]),
        int(scan_cfg["n_doses"]),
    )


def run_dose_scan(
    base_config: Dict,
    t_end: float | None = None,
    dt: float | None = None,
) -> List[Dict[str, Any]]:
    """
    Run dose scan across a configured dose grid.
    """
    dose_grid = build_dose_grid(base_config)

    if t_end is None:
        t_end = float(base_config["simulation"]["default_t_end"])
    if dt is None:
        dt = float(base_config["simulation"]["default_dt"])

    scan_results = []

    for dose in dose_grid:
        config = copy.deepcopy(base_config)

        sim = run_simulation(
            config=config,
            dose=float(dose),
            t_end=t_end,
            dt=dt,
        )

        summary = extract_summary_metrics(sim)
        feasibility = evaluate_therapeutic_window(summary, config)

        row = {
            "dose": float(dose),
            **summary,
            **feasibility,
        }
        scan_results.append(row)

    return scan_results


def find_min_effective_dose(scan_results: List[Dict[str, Any]]) -> Dict[str, Any] | None:
    """
    Find the minimum feasible dose from dose scan results.
    """
    feasible_rows = [row for row in scan_results if row["feasible"]]
    if not feasible_rows:
        return None
    feasible_rows = sorted(feasible_rows, key=lambda x: x["dose"])
    return feasible_rows[0]