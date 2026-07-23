from __future__ import annotations

from typing import Dict
import numpy as np


def validate_editing_params(params: Dict[str, float]) -> None:
    required = [
        "k_cat_on",
        "K_d_on",
        "k_cat_off",
        "K_d_off",
        "k_loss_on",
        "k_loss_off",
    ]
    for key in required:
        if key not in params:
            raise KeyError(f"Missing editing parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"Editing parameter '{key}' must be non-negative.")


def compute_editing_fluxes(
    P_brain: float,
    S_on: float,
    S_off: float,
    E_on: float,
    E_off: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Competitive editing fluxes for module 5.
    """
    k_cat_on = params["k_cat_on"]
    K_d_on = params["K_d_on"]
    k_cat_off = params["k_cat_off"]
    K_d_off = params["K_d_off"]
    k_loss_on = params["k_loss_on"]
    k_loss_off = params["k_loss_off"]

    # Competitive catalytic forms
    v_on = (
        k_cat_on * P_brain * S_on
        / (K_d_on * (1.0 + S_off / max(K_d_off, 1e-12)) + S_on + 1e-12)
    )

    v_off = (
        k_cat_off * P_brain * S_off
        / (K_d_off * (1.0 + S_on / max(K_d_on, 1e-12)) + S_off + 1e-12)
    )

    fluxes = {
        "v_on": v_on,
        "v_off": v_off,
        "loss_on": k_loss_on * E_on,
        "loss_off": k_loss_off * E_off,
    }
    return fluxes


def competitive_editing_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Module 5: competitive on-target / off-target editing dynamics.

    States affected:
    - S_on
    - S_off
    - E_on
    - E_off
    """
    _ = t
    validate_editing_params(params)

    dydt = np.zeros_like(y)

    P_brain = y[idx["P_brain"]]
    S_on = y[idx["S_on"]]
    S_off = y[idx["S_off"]]
    E_on = y[idx["E_on"]]
    E_off = y[idx["E_off"]]

    fluxes = compute_editing_fluxes(
        P_brain=P_brain,
        S_on=S_on,
        S_off=S_off,
        E_on=E_on,
        E_off=E_off,
        params=params,
    )

    v_on = fluxes["v_on"]
    v_off = fluxes["v_off"]
    loss_on = fluxes["loss_on"]
    loss_off = fluxes["loss_off"]

    # substrate depletion
    dydt[idx["S_on"]] += -v_on
    dydt[idx["S_off"]] += -v_off

    # edited outcome accumulation
    dydt[idx["E_on"]] += v_on - loss_on
    dydt[idx["E_off"]] += v_off - loss_off

    return dydt


def extract_editing_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct module-5 editing flux trajectories from simulation result.
    """
    flux_history = {
        "v_on": [],
        "v_off": [],
        "loss_on": [],
        "loss_off": [],
    }

    for k in range(y_array.shape[1]):
        P_brain = y_array[idx["P_brain"], k]
        S_on = y_array[idx["S_on"], k]
        S_off = y_array[idx["S_off"], k]
        E_on = y_array[idx["E_on"], k]
        E_off = y_array[idx["E_off"], k]

        fluxes = compute_editing_fluxes(
            P_brain=P_brain,
            S_on=S_on,
            S_off=S_off,
            E_on=E_on,
            E_off=E_off,
            params=params,
        )

        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history


def compute_editing_metrics(
    y: np.ndarray,
    idx: Dict[str, int],
    config: Dict,
) -> Dict[str, np.ndarray]:
    """
    Compute derived module-5 metrics:
    - on-target editing rate
    - off-target burden
    - specificity index
    """
    S_on_init = float(config["editing"]["S_on_init"])
    S_off_init = float(config["editing"]["S_off_init"])
    eps = float(config["editing"].get("specificity_eps", 1e-12))

    E_on = y[idx["E_on"], :]
    E_off = y[idx["E_off"], :]

    on_target_editing_rate = E_on / max(S_on_init, eps)
    off_target_burden = E_off / max(S_off_init, eps)
    specificity_index = E_on / (E_off + eps)

    return {
        "on_target_editing_rate": on_target_editing_rate,
        "off_target_burden": off_target_burden,
        "specificity_index": specificity_index,
    }