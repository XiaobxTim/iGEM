from __future__ import annotations

from typing import Dict
import numpy as np


REQUIRED_KEYS = (
    "editor_type",
    "k_on_112",
    "k_off_112",
    "k_cat_112",
    "k_on_158",
    "k_off_158",
    "k_cat_158",
    "k_prod_apoe",
    "k_deg_apoe",
    "local_bystander_per_112",
    "local_bystander_per_158",
    "k_on_puf_off",
    "k_off_puf_off",
    "k_cat_puf_off",
    "k_deaminase_bg",
    "k_prod_puf_off",
    "k_deg_puf_off",
    "k_prod_deaminase_bg",
    "k_deg_deaminase_bg",
)


def validate_params(params: Dict) -> None:
    for key in REQUIRED_KEYS:
        if key not in params:
            raise KeyError(f"Missing APOE multisite editing parameter: {key}")
        if key != "editor_type" and float(params[key]) < 0:
            raise ValueError(f"APOE multisite parameter '{key}' must be non-negative.")


def _positive(value: float) -> float:
    return max(float(value), 0.0)


def _scale_for_editor_type(params: Dict) -> Dict[str, float]:
    """
    Return editor-specific modifiers.

    v4 explicitly distinguishes APOBEC3A/A3A-like editors from APOBEC1-like
    editors. These are not final literature-calibrated constants; they are
    transparent priors that prevent us from mixing A3A and APOBEC1 mechanisms.
    """
    editor_type = str(params.get("editor_type", "A3A")).lower()
    if editor_type in {"a3a", "apobec3a", "engineered_a3a"}:
        return {
            "cat_scale": 1.0,
            "uc_context_scale": float(params.get("uc_context_scale", 1.2)),
            "background_scale": float(params.get("a3a_background_scale", 1.0)),
        }
    if editor_type in {"apobec1", "apobec1_like"}:
        return {
            "cat_scale": float(params.get("apobec1_cat_scale", 0.7)),
            "uc_context_scale": float(params.get("apobec1_uc_context_scale", 0.9)),
            "background_scale": float(params.get("apobec1_background_scale", 0.4)),
        }
    if editor_type in {"proapobec", "proapobec_like"}:
        return {
            "cat_scale": float(params.get("proapobec_cat_scale", 1.25)),
            "uc_context_scale": float(params.get("proapobec_uc_context_scale", 1.1)),
            "background_scale": float(params.get("proapobec_background_scale", 0.5)),
        }
    raise ValueError(f"Unsupported editor_type '{params.get('editor_type')}'.")


def compute_apoe_multisite_fluxes(y: np.ndarray, params: Dict, idx: Dict[str, int]) -> Dict[str, float]:
    """
    Compute v4 APOE112/APOE158 and three-class off-target fluxes.

    APOE biology encoded here:
    - APOE4 has Arg112/Arg158.
    - Editing 112 only gives APOE3-like.
    - Editing both 112 and 158 gives APOE2-like.
    - Editing 158 only is tracked as an incomplete/mixed edit because it is not
      the normal APOE3-like therapeutic objective.
    """
    validate_params(params)
    editor_mod = _scale_for_editor_type(params)

    p = _positive(y[idx["P_brain"]])
    s4 = _positive(y[idx["S_APOE4"]])
    s3 = _positive(y[idx["S_APOE3_like"]])
    s2 = _positive(y[idx["S_APOE2_like"]])
    s158_only = _positive(y[idx["S_APOE158_only"]])
    c112 = _positive(y[idx["C_APOE112"]])
    c158 = _positive(y[idx["C_APOE158"]])
    c158_after112 = _positive(y[idx["C_APOE158_after112"]])
    c112_after158 = _positive(y[idx["C_APOE112_after158"]])
    s_puf_off = _positive(y[idx["S_puf_off"]])
    c_puf_off = _positive(y[idx["C_puf_off"]])
    s_deam = _positive(y[idx["S_deaminase_bg"]])

    cat112 = params["k_cat_112"] * editor_mod["cat_scale"] * editor_mod["uc_context_scale"]
    cat158 = params["k_cat_158"] * editor_mod["cat_scale"]
    cat_off = params["k_cat_puf_off"] * editor_mod["cat_scale"]
    k_bg = params["k_deaminase_bg"] * editor_mod["background_scale"]

    bind112 = params["k_on_112"] * p * s4
    unbind112 = params["k_off_112"] * c112
    edit112 = cat112 * c112

    bind158 = params["k_on_158"] * p * s4
    unbind158 = params["k_off_158"] * c158
    edit158 = cat158 * c158

    bind158_after112 = params["k_on_158"] * p * s3
    unbind158_after112 = params["k_off_158"] * c158_after112
    edit158_after112 = cat158 * c158_after112

    bind112_after158 = params["k_on_112"] * p * s158_only
    unbind112_after158 = params["k_off_112"] * c112_after158
    edit112_after158 = cat112 * c112_after158

    bind_puf_off = params["k_on_puf_off"] * p * s_puf_off
    unbind_puf_off = params["k_off_puf_off"] * c_puf_off
    edit_puf_off = cat_off * c_puf_off

    deaminase_bg_edit = k_bg * p * s_deam
    local_bystander = (
        params["local_bystander_per_112"] * (edit112 + edit112_after158)
        + params["local_bystander_per_158"] * (edit158 + edit158_after112)
    )

    return {
        "prod_apoe": params["k_prod_apoe"],
        "deg_apoe4": params["k_deg_apoe"] * s4,
        "deg_apoe3_like": params["k_deg_apoe"] * s3,
        "deg_apoe2_like": params["k_deg_apoe"] * s2,
        "deg_apoe158_only": params["k_deg_apoe"] * s158_only,
        "bind112": bind112,
        "unbind112": unbind112,
        "edit112": edit112,
        "bind158": bind158,
        "unbind158": unbind158,
        "edit158": edit158,
        "bind158_after112": bind158_after112,
        "unbind158_after112": unbind158_after112,
        "edit158_after112": edit158_after112,
        "bind112_after158": bind112_after158,
        "unbind112_after158": unbind112_after158,
        "edit112_after158": edit112_after158,
        "local_bystander": local_bystander,
        "prod_puf_off": params["k_prod_puf_off"],
        "deg_puf_off": params["k_deg_puf_off"] * s_puf_off,
        "bind_puf_off": bind_puf_off,
        "unbind_puf_off": unbind_puf_off,
        "edit_puf_off": edit_puf_off,
        "prod_deaminase_bg": params["k_prod_deaminase_bg"],
        "deg_deaminase_bg": params["k_deg_deaminase_bg"] * s_deam,
        "edit_deaminase_bg": deaminase_bg_edit,
    }


def apoe_multisite_rhs(t: float, y: np.ndarray, params: Dict, idx: Dict[str, int]) -> np.ndarray:
    _ = t
    fluxes = compute_apoe_multisite_fluxes(y, params, idx)
    dydt = np.zeros_like(y)

    editor_binding_loss = (
        fluxes["bind112"] + fluxes["bind158"]
        + fluxes["bind158_after112"] + fluxes["bind112_after158"]
        + fluxes["bind_puf_off"]
    )
    editor_release = (
        fluxes["unbind112"] + fluxes["unbind158"]
        + fluxes["unbind158_after112"] + fluxes["unbind112_after158"]
        + fluxes["unbind_puf_off"]
        + fluxes["edit112"] + fluxes["edit158"]
        + fluxes["edit158_after112"] + fluxes["edit112_after158"]
        + fluxes["edit_puf_off"]
    )
    dydt[idx["P_brain"]] += -editor_binding_loss + editor_release

    dydt[idx["S_APOE4"]] += (
        fluxes["prod_apoe"] - fluxes["deg_apoe4"]
        - fluxes["bind112"] + fluxes["unbind112"]
        - fluxes["bind158"] + fluxes["unbind158"]
    )
    dydt[idx["C_APOE112"]] += fluxes["bind112"] - fluxes["unbind112"] - fluxes["edit112"]
    dydt[idx["C_APOE158"]] += fluxes["bind158"] - fluxes["unbind158"] - fluxes["edit158"]

    dydt[idx["S_APOE3_like"]] += (
        fluxes["edit112"] - fluxes["deg_apoe3_like"]
        - fluxes["bind158_after112"] + fluxes["unbind158_after112"]
    )
    dydt[idx["C_APOE158_after112"]] += (
        fluxes["bind158_after112"] - fluxes["unbind158_after112"] - fluxes["edit158_after112"]
    )

    dydt[idx["S_APOE158_only"]] += (
        fluxes["edit158"] - fluxes["deg_apoe158_only"]
        - fluxes["bind112_after158"] + fluxes["unbind112_after158"]
    )
    dydt[idx["C_APOE112_after158"]] += (
        fluxes["bind112_after158"] - fluxes["unbind112_after158"] - fluxes["edit112_after158"]
    )

    dydt[idx["S_APOE2_like"]] += (
        fluxes["edit158_after112"] + fluxes["edit112_after158"] - fluxes["deg_apoe2_like"]
    )

    dydt[idx["B_local_bystander"]] += fluxes["local_bystander"]

    dydt[idx["S_puf_off"]] += (
        fluxes["prod_puf_off"] - fluxes["deg_puf_off"]
        - fluxes["bind_puf_off"] + fluxes["unbind_puf_off"]
    )
    dydt[idx["C_puf_off"]] += (
        fluxes["bind_puf_off"] - fluxes["unbind_puf_off"] - fluxes["edit_puf_off"]
    )
    dydt[idx["E_puf_off"]] += fluxes["edit_puf_off"]

    dydt[idx["S_deaminase_bg"]] += (
        fluxes["prod_deaminase_bg"]
        - fluxes["deg_deaminase_bg"]
        - fluxes["edit_deaminase_bg"]
    )
    dydt[idx["E_deaminase_bg"]] += fluxes["edit_deaminase_bg"]

    # Compatibility accumulators used by older optimization/reporting code.
    dydt[idx["E_on"]] += fluxes["edit112"] + fluxes["edit158_after112"]
    dydt[idx["E_off"]] += (
        fluxes["local_bystander"]
        + fluxes["edit_puf_off"]
        + fluxes["edit_deaminase_bg"]
    )

    return dydt


def extract_apoe_multisite_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict,
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    _ = t_array
    history: Dict[str, list[float]] = {}
    for col in range(y_array.shape[1]):
        fluxes = compute_apoe_multisite_fluxes(y_array[:, col], params, idx)
        for key, value in fluxes.items():
            history.setdefault(key, []).append(value)
    return {key: np.asarray(values, dtype=float) for key, values in history.items()}


def compute_apoe_multisite_metrics(y: np.ndarray, idx: Dict[str, int], config: Dict) -> Dict[str, np.ndarray]:
    eps = float(config["editing"].get("specificity_eps", 1e-12))
    apoe4 = y[idx["S_APOE4"], :]
    apoe3 = y[idx["S_APOE3_like"], :]
    apoe2 = y[idx["S_APOE2_like"], :]
    apoe158_only = y[idx["S_APOE158_only"], :]
    complexes = (
        y[idx["C_APOE112"], :]
        + y[idx["C_APOE158"], :]
        + y[idx["C_APOE158_after112"], :]
        + y[idx["C_APOE112_after158"], :]
    )
    total_apoe = apoe4 + apoe3 + apoe2 + apoe158_only + complexes + eps

    local = y[idx["B_local_bystander"], :]
    puf_off = y[idx["E_puf_off"], :]
    deam = y[idx["E_deaminase_bg"], :]
    off_total = local + puf_off + deam

    therapeutic = apoe3 + apoe2
    apoe3_like_fraction = apoe3 / total_apoe
    apoe2_like_fraction = apoe2 / total_apoe
    apoe4_fraction = apoe4 / total_apoe
    mixed_edit_fraction = apoe158_only / total_apoe

    return {
        "on_target_editing_rate": therapeutic / total_apoe,
        "off_target_burden": off_total / total_apoe,
        "specificity_index": therapeutic / (off_total + eps),
        "on_target_editing_fraction": therapeutic / total_apoe,
        "off_target_editing_fraction": off_total / (total_apoe + off_total + eps),
        "apoe4_fraction": apoe4_fraction,
        "apoe3_like_fraction": apoe3_like_fraction,
        "apoe2_like_fraction": apoe2_like_fraction,
        "mixed_edit_fraction": mixed_edit_fraction,
        "local_bystander_burden": local / total_apoe,
        "puf_mismatch_burden": puf_off / total_apoe,
        "deaminase_background_burden": deam / total_apoe,
        "ldlr_binding_risk_proxy": apoe2 / total_apoe,
    }
