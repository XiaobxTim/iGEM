from __future__ import annotations

from typing import Dict


def evaluate_therapeutic_window(
    summary_metrics: Dict[str, float],
    config: Dict,
) -> Dict[str, object]:
    """
    Determine whether a simulation result is feasible under current thresholds.
    """
    opt_cfg = config["optimization"]

    on_target_ok = (
        summary_metrics["on_target_editing_rate_final"]
        >= opt_cfg["E_on_target_threshold"]
    )

    off_target_ok = (
        summary_metrics["off_target_burden_final"]
        <= opt_cfg["E_off_max"]
    )

    specificity_ok = (
        summary_metrics["specificity_index_final"]
        >= opt_cfg["SI_min"]
    )

    liver_ok = (
        summary_metrics["AUC_liver"]
        <= opt_cfg["AUC_liver_max"]
    )

    blood_ok = (
        summary_metrics["Cmax_blood"]
        <= opt_cfg["Cmax_blood_max"]
    )

    feasible = all([
        on_target_ok,
        off_target_ok,
        specificity_ok,
        liver_ok,
        blood_ok,
    ])

    return {
        "feasible": feasible,
        "on_target_ok": on_target_ok,
        "off_target_ok": off_target_ok,
        "specificity_ok": specificity_ok,
        "liver_ok": liver_ok,
        "blood_ok": blood_ok,
    }