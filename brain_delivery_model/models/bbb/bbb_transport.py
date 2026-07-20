from __future__ import annotations

from typing import Dict
import numpy as np


def validate_bbb_params(params: Dict[str, float]) -> None:
    required = [
        "k_brainblood_to_EC",
        "k_EC_to_endo",
        "k_EC_to_brainblood",
        "k_endo_to_ISF",
        "k_endo_to_brainblood",
        "k_endo_loss",
        "k_clear_ISF",
    ]
    for key in required:
        if key not in params:
            raise KeyError(f"Missing BBB parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"BBB parameter '{key}' must be non-negative.")


def compute_bbb_fluxes(
    A_brain_blood: float,
    A_brain_EC: float,
    A_brain_endo: float,
    A_brain_ISF: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute all module-3 BBB transport fluxes.
    """
    fluxes = {
        "brainblood_to_EC": params["k_brainblood_to_EC"] * A_brain_blood,
        "EC_to_endo": params["k_EC_to_endo"] * A_brain_EC,
        "EC_to_brainblood": params["k_EC_to_brainblood"] * A_brain_EC,

        "endo_to_ISF": params["k_endo_to_ISF"] * A_brain_endo,
        "endo_to_brainblood": params["k_endo_to_brainblood"] * A_brain_endo,
        "endo_loss": params["k_endo_loss"] * A_brain_endo,

        "clear_ISF": params["k_clear_ISF"] * A_brain_ISF,
    }
    return fluxes


def bbb_transport_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Module 3: BBB transport dynamics.

    States affected:
    - A_brain_blood
    - A_brain_EC
    - A_brain_endo
    - A_brain_ISF
    - A_cleared (for endo_loss + ISF clearance accumulation)
    """
    _ = t
    validate_bbb_params(params)

    dydt = np.zeros_like(y)

    A_brain_blood = y[idx["A_brain_blood"]]
    A_brain_EC = y[idx["A_brain_EC"]]
    A_brain_endo = y[idx["A_brain_endo"]]
    A_brain_ISF = y[idx["A_brain_ISF"]]

    fluxes = compute_bbb_fluxes(
        A_brain_blood=A_brain_blood,
        A_brain_EC=A_brain_EC,
        A_brain_endo=A_brain_endo,
        A_brain_ISF=A_brain_ISF,
        params=params,
    )

    # A_brain_blood:
    # lose to EC, gain from EC dissociation and endo recycling
    dydt[idx["A_brain_blood"]] += (
        - fluxes["brainblood_to_EC"]
        + fluxes["EC_to_brainblood"]
        + fluxes["endo_to_brainblood"]
    )

    # A_brain_EC:
    # gain from brain blood, lose to endo and back to brain blood
    dydt[idx["A_brain_EC"]] += (
        fluxes["brainblood_to_EC"]
        - fluxes["EC_to_endo"]
        - fluxes["EC_to_brainblood"]
    )

    # A_brain_endo:
    # gain from EC, lose to ISF, brainblood recycling, and degradation
    dydt[idx["A_brain_endo"]] += (
        fluxes["EC_to_endo"]
        - fluxes["endo_to_ISF"]
        - fluxes["endo_to_brainblood"]
        - fluxes["endo_loss"]
    )

    # A_brain_ISF:
    # gain from successful transcytosis, lose by local clearance
    # later module 4 will additionally remove from ISF via cell uptake
    dydt[idx["A_brain_ISF"]] += (
        fluxes["endo_to_ISF"]
        - fluxes["clear_ISF"]
    )

    # Cleared accumulator
    dydt[idx["A_cleared"]] += (
        fluxes["endo_loss"]
        + fluxes["clear_ISF"]
    )

    return dydt


def extract_bbb_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct module-3 BBB flux trajectories from simulation result.
    """
    flux_history = {
        "brainblood_to_EC": [],
        "EC_to_endo": [],
        "EC_to_brainblood": [],
        "endo_to_ISF": [],
        "endo_to_brainblood": [],
        "endo_loss": [],
        "clear_ISF": [],
    }

    for k in range(y_array.shape[1]):
        A_brain_blood = y_array[idx["A_brain_blood"], k]
        A_brain_EC = y_array[idx["A_brain_EC"], k]
        A_brain_endo = y_array[idx["A_brain_endo"], k]
        A_brain_ISF = y_array[idx["A_brain_ISF"], k]

        fluxes = compute_bbb_fluxes(
            A_brain_blood=A_brain_blood,
            A_brain_EC=A_brain_EC,
            A_brain_endo=A_brain_endo,
            A_brain_ISF=A_brain_ISF,
            params=params,
        )

        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history