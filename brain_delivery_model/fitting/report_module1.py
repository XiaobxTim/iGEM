from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import json

import pandas as pd

from fitting.objective_module1 import SOFT_RANGES
from fitting.plot_module1_fit import (
    plot_module1_fit_observations,
    plot_module1_fit_states,
    plot_module1_fit_fluxes,
)


def _ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def build_soft_range_report(best_params: Dict[str, float]) -> pd.DataFrame:
    """
    Build a dataframe showing whether each fitted parameter falls within soft range.
    """
    rows = []
    for name, value in best_params.items():
        low, high = SOFT_RANGES[name]
        in_range = (low <= value <= high)
        rows.append({
            "parameter": name,
            "value": value,
            "soft_low": low,
            "soft_high": high,
            "in_soft_range": in_range,
        })
    return pd.DataFrame(rows)


def build_parameter_table(best_params: Dict[str, float]) -> pd.DataFrame:
    """
    Convert fitted parameter dictionary to dataframe.
    """
    return pd.DataFrame(
        [{"parameter": k, "value": v} for k, v in best_params.items()]
    )


def save_module1_report(
    fit_result: Dict[str, Any],
    dose: float,
    observations: Dict,
    output_dir: str | Path,
    t_end: float = 72.0,
) -> Dict[str, Path]:
    """
    Save module 1 fitting report including:
    - parameter table CSV
    - soft range check CSV
    - summary JSON
    - fit figures

    Parameters
    ----------
    fit_result : dict
        Output of fit_module1_two_stage(...)
    dose : float
        Initial dose.
    observations : dict
        Observation dictionary.
    output_dir : str or Path
        Output directory for report.
    t_end : float
        Plotting horizon.

    Returns
    -------
    dict
        Paths of saved outputs.
    """
    output_dir = _ensure_dir(output_dir)
    best_params = fit_result["best_params"]

    # Tables
    param_df = build_parameter_table(best_params)
    soft_df = build_soft_range_report(best_params)

    param_csv = output_dir / "module1_best_params.csv"
    soft_csv = output_dir / "module1_soft_range_check.csv"

    param_df.to_csv(param_csv, index=False)
    soft_df.to_csv(soft_csv, index=False)

    # Summary JSON
    summary = {
        "best_params": best_params,
        "best_loss": fit_result["best_loss"],
        "global_result": {
            "loss": fit_result["global_result"]["loss"],
            "success": fit_result["global_result"]["success"],
            "message": str(fit_result["global_result"]["message"]),
            "method": fit_result["global_result"]["method"],
        },
        "local_result": {
            "loss": fit_result["local_result"]["loss"],
            "success": fit_result["local_result"]["success"],
            "message": str(fit_result["local_result"]["message"]),
            "method": fit_result["local_result"]["method"],
        },
    }

    summary_json = output_dir / "module1_fit_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # Figures
    fit_obs_png = plot_module1_fit_observations(
        best_params=best_params,
        dose=dose,
        observations=observations,
        output_dir=output_dir,
        t_end=t_end,
        filename="module1_fit_observations.png",
    )

    fit_states_png = plot_module1_fit_states(
        best_params=best_params,
        dose=dose,
        output_dir=output_dir,
        t_end=t_end,
        filename="module1_fit_states.png",
    )

    fit_fluxes_png = plot_module1_fit_fluxes(
        best_params=best_params,
        dose=dose,
        output_dir=output_dir,
        t_end=t_end,
        filename="module1_fit_fluxes.png",
    )

    return {
        "param_csv": param_csv,
        "soft_csv": soft_csv,
        "summary_json": summary_json,
        "fit_obs_png": fit_obs_png,
        "fit_states_png": fit_states_png,
        "fit_fluxes_png": fit_fluxes_png,
    }