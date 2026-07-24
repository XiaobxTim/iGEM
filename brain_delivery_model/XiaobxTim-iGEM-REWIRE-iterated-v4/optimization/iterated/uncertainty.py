from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, Tuple
import csv
import math
import random

import yaml

from models.full_model.simulator import run_simulation
from optimization.metrics import extract_summary_metrics


def _get_nested(config: Dict, path: str):
    node = config
    parts = path.split(".")
    for part in parts:
        node = node[part]
    return node


def _set_nested(config: Dict, path: str, value: float) -> None:
    node = config
    parts = path.split(".")
    for part in parts[:-1]:
        node = node[part]
    node[parts[-1]] = float(value)


def load_parameter_priors(path: str | Path) -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data.get("parameters", {})


def sample_parameter(prior: Dict, rng: random.Random) -> float:
    low = float(prior["low"])
    high = float(prior["high"])
    distribution = prior.get("distribution", "uniform")
    if distribution == "log_uniform":
        return math.exp(rng.uniform(math.log(low), math.log(high)))
    if distribution == "uniform":
        return rng.uniform(low, high)
    raise ValueError(f"Unsupported prior distribution: {distribution}")


def sample_config(base_config: Dict, priors: Dict, rng: random.Random) -> Tuple[Dict, Dict[str, float]]:
    """
    Draw one uncertain parameter set and apply it to a copy of the config.
    """
    config = deepcopy(base_config)
    sampled = {}
    for path, prior in priors.items():
        if prior.get("enabled", True) is False:
            continue
        # Skip missing paths gracefully so the provenance table can include
        # future parameters before the model code uses them.
        try:
            _get_nested(config, path)
        except KeyError:
            continue
        value = sample_parameter(prior, rng)
        _set_nested(config, path, value)
        sampled[path] = value
    return config, sampled


def run_uncertainty_scan(
    base_config: Dict,
    priors: Dict,
    n_samples: int,
    dose: float,
    t_end: float,
    dt: float,
    seed: int = 7,
) -> Iterable[Dict[str, float]]:
    """
    Monte Carlo scan over uncertain biological parameters.

    The output rows include both sampled parameters and decision metrics. This
    makes the analysis useful for iGEM-style "what should we measure next?"
    reasoning, because high-spread metrics point to weakly constrained biology.
    """
    rng = random.Random(seed)
    for sample_id in range(int(n_samples)):
        config, sampled = sample_config(base_config, priors, rng)
        sim = run_simulation(config=config, dose=dose, t_end=t_end, dt=dt)
        summary = extract_summary_metrics(sim)
        yield {
            "sample_id": sample_id,
            "dose": dose,
            **sampled,
            **summary,
        }


def save_rows(rows: Iterable[Dict], output_path: str | Path) -> Path:
    rows = list(rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path
