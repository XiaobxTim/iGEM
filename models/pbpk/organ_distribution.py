from __future__ import annotations

from typing import Dict
import numpy as np


def validate_distribution_params(params: Dict[str, float]) -> None:
    required = [
        "k_blood_to_liver",
        "k_blood_to_peripheral",
        "k_blood_to_brain",
        "k_liver_to_blood",
        "k_peripheral_to_blood",
        "k_brain_to_blood",
        "k_clear_blood",
        "k_clear_liver",
        "k_clear_peripheral",
        "k_clear_brainblood",
    ]
    for key in required:
        if key not in params:
            raise KeyError(f"Missing distribution parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"Distribution parameter '{key}' must be non-negative.")


def compute_distribution_fluxes(
    A_blood: float,
    A_liver: float,
    A_peripheral: float,
    A_brain_blood: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute all module-2 fluxes.
    """
    fluxes = {
        "blood_to_liver": params["k_blood_to_liver"] * A_blood,
        "blood_to_peripheral": params["k_blood_to_peripheral"] * A_blood,
        "blood_to_brain": params["k_blood_to_brain"] * A_blood,

        "liver_to_blood": params["k_liver_to_blood"] * A_liver,
        "peripheral_to_blood": params["k_peripheral_to_blood"] * A_peripheral,
        "brain_to_blood": params["k_brain_to_blood"] * A_brain_blood,

        "clear_blood": params["k_clear_blood"] * A_blood,
        "clear_liver": params["k_clear_liver"] * A_liver,
        "clear_peripheral": params["k_clear_peripheral"] * A_peripheral,
        "clear_brainblood": params["k_clear_brainblood"] * A_brain_blood,
    }
    return fluxes


def organ_distribution_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Module 2: whole-body distribution
    blood <-> liver
    blood <-> peripheral
    blood <-> brain_blood
    + clearance
    """
    _ = t
    validate_distribution_params(params)

    dydt = np.zeros_like(y)

    A_blood = y[idx["A_blood"]]
    A_liver = y[idx["A_liver"]]
    A_peripheral = y[idx["A_peripheral"]]
    A_brain_blood = y[idx["A_brain_blood"]]

    fluxes = compute_distribution_fluxes(
        A_blood=A_blood,
        A_liver=A_liver,
        A_peripheral=A_peripheral,
        A_brain_blood=A_brain_blood,
        params=params,
    )

    # Blood
    dydt[idx["A_blood"]] += (
        - fluxes["blood_to_liver"]
        - fluxes["blood_to_peripheral"]
        - fluxes["blood_to_brain"]
        - fluxes["clear_blood"]
        + fluxes["liver_to_blood"]
        + fluxes["peripheral_to_blood"]
        + fluxes["brain_to_blood"]
    )

    # Liver
    dydt[idx["A_liver"]] += (
        fluxes["blood_to_liver"]
        - fluxes["liver_to_blood"]
        - fluxes["clear_liver"]
    )

    # Peripheral lumped compartment
    dydt[idx["A_peripheral"]] += (
        fluxes["blood_to_peripheral"]
        - fluxes["peripheral_to_blood"]
        - fluxes["clear_peripheral"]
    )

    # Brain vascular side
    dydt[idx["A_brain_blood"]] += (
        fluxes["blood_to_brain"]
        - fluxes["brain_to_blood"]
        - fluxes["clear_brainblood"]
    )

    # Cleared accumulator
    dydt[idx["A_cleared"]] += (
        fluxes["clear_blood"]
        + fluxes["clear_liver"]
        + fluxes["clear_peripheral"]
        + fluxes["clear_brainblood"]
    )

    return dydt


def extract_distribution_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct module-2 flux trajectories from simulation result.
    """
    flux_history = {
        "blood_to_liver": [],
        "blood_to_peripheral": [],
        "blood_to_brain": [],
        "liver_to_blood": [],
        "peripheral_to_blood": [],
        "brain_to_blood": [],
        "clear_blood": [],
        "clear_liver": [],
        "clear_peripheral": [],
        "clear_brainblood": [],
    }

    for k in range(y_array.shape[1]):
        A_blood = y_array[idx["A_blood"], k]
        A_liver = y_array[idx["A_liver"], k]
        A_peripheral = y_array[idx["A_peripheral"], k]
        A_brain_blood = y_array[idx["A_brain_blood"], k]

        fluxes = compute_distribution_fluxes(
            A_blood=A_blood,
            A_liver=A_liver,
            A_peripheral=A_peripheral,
            A_brain_blood=A_brain_blood,
            params=params,
        )

        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history