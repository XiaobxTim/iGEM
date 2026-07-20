from __future__ import annotations

from typing import Dict
import numpy as np


REQUIRED_ABSORPTION_KEYS = (
    "k_lymph",
    "k_blood",
    "k_deg_loc",
    "k_lymph_to_blood",
)


def get_absorption_params(config: Dict, route: str) -> Dict[str, float]:
    """
    Extract route-specific absorption parameters from the global configuration.
    """
    if "absorption" not in config:
        raise KeyError("Missing 'absorption' section in config.")
    if route not in config["absorption"]:
        raise KeyError(f"Unsupported route '{route}' in absorption config.")
    params = config["absorption"][route]
    validate_absorption_params(params)
    return params


def validate_absorption_params(params: Dict[str, float]) -> None:
    """Validate that all absorption parameters exist and are non-negative."""
    for key in REQUIRED_ABSORPTION_KEYS:
        if key not in params:
            raise KeyError(f"Missing absorption parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"Absorption parameter '{key}' must be non-negative.")


def compute_absorption_fluxes(
    A_dep: float,
    A_lymph: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute the four physical fluxes in module 1.

    Returns
    -------
    dict
        Keys:
        - to_lymph
        - to_blood
        - local_loss
        - lymph_to_blood
    """
    validate_absorption_params(params)

    k_lymph = params["k_lymph"]
    k_blood = params["k_blood"]
    k_deg_loc = params["k_deg_loc"]
    k_lymph_to_blood = params["k_lymph_to_blood"]

    return {
        "to_lymph": k_lymph * A_dep,
        "to_blood": k_blood * A_dep,
        "local_loss": k_deg_loc * A_dep,
        "lymph_to_blood": k_lymph_to_blood * A_lymph,
    }


def lymphatic_absorption_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    RHS contribution of the absorption module.

    This module governs:
        depot -> lymph
        depot -> blood
        depot -> local loss
        lymph -> blood
    """
    _ = t  # reserved for future time-dependent extensions

    dydt = np.zeros_like(y)

    A_dep = y[idx["A_dep"]]
    A_lymph = y[idx["A_lymph"]]

    fluxes = compute_absorption_fluxes(A_dep=A_dep, A_lymph=A_lymph, params=params)

    to_lymph = fluxes["to_lymph"]
    to_blood = fluxes["to_blood"]
    local_loss = fluxes["local_loss"]
    lymph_to_blood = fluxes["lymph_to_blood"]

    dydt[idx["A_dep"]] += -(to_lymph + to_blood + local_loss)
    dydt[idx["A_lymph"]] += to_lymph - lymph_to_blood
    dydt[idx["A_blood"]] += to_blood + lymph_to_blood

    return dydt


def extract_absorption_flux_trajectories(
    x_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct module-1 fluxes from a simulated trajectory.

    Parameters
    ----------
    y_matrix : np.ndarray
        Shape (n_states, n_times).
    params : dict
        Absorption parameter dictionary.
    idx : dict
        State index map.
    """
    flux_history = {
        "to_lymph": [],
        "to_blood": [],
        "local_loss": [],
        "lymph_to_blood": [],
    }

    for k in range(y_array.shape[1]):
        A_dep = float(y_array[idx["A_dep"], k])
        A_lymph = float(y_array[idx["A_lymph"], k])
        fluxes = compute_absorption_fluxes(A_dep=A_dep, A_lymph=A_lymph, params=params)
        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history


__all__ = [
    "get_absorption_params",
    "validate_absorption_params",
    "compute_absorption_fluxes",
    "lymphatic_absorption_rhs",
    "extract_absorption_flux_trajectories",
]
