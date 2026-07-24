from __future__ import annotations

from typing import Dict, Any
import numpy as np

from models.editing.module5 import compute_editing_metrics


def trapz_auc(t: np.ndarray, y: np.ndarray) -> float:
    return float(np.trapz(y, t))


def extract_summary_metrics(results: Dict[str, Any]) -> Dict[str, float]:
    """
    Extract summary metrics from a full simulation result.

    Returns
    -------
    dict
        Summary metrics for module 6.
    """
    t = results["t"]
    y = results["y"]
    idx = results["idx"]
    config = results["config"]

    metrics = compute_editing_metrics(y=y, idx=idx, config=config)

    E_on_final = float(y[idx["E_on"], -1])
    E_off_final = float(y[idx["E_off"], -1])

    on_target_editing_rate_final = float(metrics["on_target_editing_rate"][-1])
    off_target_burden_final = float(metrics["off_target_burden"][-1])
    specificity_index_final = float(metrics["specificity_index"][-1])

    P_brain_peak = float(np.max(y[idx["P_brain"], :]))
    AUC_liver = float(trapz_auc(t, y[idx["A_liver"], :]))
    Cmax_blood = float(np.max(y[idx["A_blood"], :]))

    summary = {
        "E_on_final": E_on_final,
        "E_off_final": E_off_final,
        "on_target_editing_rate_final": on_target_editing_rate_final,
        "off_target_burden_final": off_target_burden_final,
        "specificity_index_final": specificity_index_final,
        "P_brain_peak": P_brain_peak,
        "AUC_liver": AUC_liver,
        "Cmax_blood": Cmax_blood,
    }
    for name, arr in metrics.items():
        if name not in {
            "on_target_editing_rate",
            "off_target_burden",
            "specificity_index",
        }:
            summary[f"{name}_final"] = float(arr[-1])
    return summary
