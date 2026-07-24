from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, List
import csv
import math
import random

import yaml

from models.full_model.simulator import run_simulation
from models.calibration.wetlab_bridge import WetlabObservation, weighted_rmse
from optimization.iterated.uncertainty import sample_parameter


DEFAULT_FIT_PARAMETERS = [
    "distribution.k_blood_to_liver",
    "distribution.k_blood_to_brain",
    "intracellular.k_tx",
    "intracellular.k_tl",
    "intracellular.k_deg_p",
    "editing.k_cat_112",
    "editing.k_cat_158",
    "editing.k_on_puf_off",
    "editing.k_cat_puf_off",
    "editing.k_deaminase_bg",
]


def _set_nested(config: Dict, path: str, value: float) -> None:
    node = config
    parts = path.split(".")
    for part in parts[:-1]:
        node = node[part]
    node[parts[-1]] = float(value)


def _get_nested(config: Dict, path: str) -> float:
    node = config
    for part in path.split("."):
        node = node[part]
    return float(node)


def evaluate_config(
    config: Dict,
    observations: Iterable[WetlabObservation],
    dose: float,
    t_end: float,
    dt: float,
) -> tuple[float, Dict]:
    sim = run_simulation(config=config, dose=dose, t_end=t_end, dt=dt)
    return weighted_rmse(sim, observations), sim


def random_search_fit(
    base_config: Dict,
    priors: Dict,
    observations: List[WetlabObservation],
    fit_parameters: List[str] | None = None,
    n_samples: int = 32,
    dose: float = 1.0,
    t_end: float = 72.0,
    dt: float = 0.2,
    seed: int = 11,
) -> tuple[Dict, List[Dict[str, float]]]:
    """
    Fit selected high-impact parameters by transparent random search.

    This is intentionally simple for iGEM use: it is easy to explain on a wiki,
    robust with sparse early wet-lab data, and can be replaced by scipy-based
    optimization once real replicate data arrive.
    """
    rng = random.Random(seed)
    fit_parameters = fit_parameters or DEFAULT_FIT_PARAMETERS
    fit_parameters = [p for p in fit_parameters if p in priors]

    baseline_score, _ = evaluate_config(base_config, observations, dose, t_end, dt)
    best_config = deepcopy(base_config)
    best_score = baseline_score
    rows: List[Dict[str, float]] = []

    baseline_row = {"sample_id": -1, "objective": baseline_score}
    for path in fit_parameters:
        baseline_row[path] = _get_nested(base_config, path)
    rows.append(baseline_row)

    for sample_id in range(int(n_samples)):
        config = deepcopy(base_config)
        row: Dict[str, float] = {"sample_id": float(sample_id)}
        for path in fit_parameters:
            value = sample_parameter(priors[path], rng)
            _set_nested(config, path, value)
            row[path] = value
        score, _ = evaluate_config(config, observations, dose, t_end, dt)
        if not math.isfinite(score):
            score = float("inf")
        row["objective"] = score
        rows.append(row)
        if score < best_score:
            best_score = score
            best_config = config

    return best_config, rows


def save_fit_results(rows: List[Dict[str, float]], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_fit_config(config: Dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return path
