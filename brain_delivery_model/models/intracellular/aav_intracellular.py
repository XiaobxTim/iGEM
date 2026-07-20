from __future__ import annotations

from typing import Dict
import numpy as np


def validate_intracellular_params(params: Dict[str, float]) -> None:
    required = [
        "k_ISF_to_cell",
        "k_cell_to_nuc",
        "k_cell_loss",
        "k_deg_v",
        "k_tx",
        "k_deg_m",
        "k_tl",
        "k_deg_p",
    ]
    for key in required:
        if key not in params:
            raise KeyError(f"Missing intracellular parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"Intracellular parameter '{key}' must be non-negative.")


def compute_intracellular_fluxes(
    A_brain_ISF: float,
    A_brain_cell: float,
    A_brain_nuc: float,
    mRNA_brain: float,
    P_brain: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute all module-4 intracellular fluxes.
    """
    fluxes = {
        "ISF_to_cell": params["k_ISF_to_cell"] * A_brain_ISF,
        "cell_to_nuc": params["k_cell_to_nuc"] * A_brain_cell,
        "cell_loss": params["k_cell_loss"] * A_brain_cell,

        "deg_v": params["k_deg_v"] * A_brain_nuc,
        "tx": params["k_tx"] * A_brain_nuc,

        "deg_m": params["k_deg_m"] * mRNA_brain,
        "tl": params["k_tl"] * mRNA_brain,

        "deg_p": params["k_deg_p"] * P_brain,
    }
    return fluxes


def aav_intracellular_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    Module 4: brain cell uptake and expression.

    States affected:
    - A_brain_ISF
    - A_brain_cell
    - A_brain_nuc
    - mRNA_brain
    - P_brain
    - A_cleared (for cell_loss + deg_v + deg_m + deg_p if you want to accumulate losses)
    """
    _ = t
    validate_intracellular_params(params)

    dydt = np.zeros_like(y)

    A_brain_ISF = y[idx["A_brain_ISF"]]
    A_brain_cell = y[idx["A_brain_cell"]]
    A_brain_nuc = y[idx["A_brain_nuc"]]
    mRNA_brain = y[idx["mRNA_brain"]]
    P_brain = y[idx["P_brain"]]

    fluxes = compute_intracellular_fluxes(
        A_brain_ISF=A_brain_ISF,
        A_brain_cell=A_brain_cell,
        A_brain_nuc=A_brain_nuc,
        mRNA_brain=mRNA_brain,
        P_brain=P_brain,
        params=params,
    )

    # A_brain_ISF:
    # lose by cell uptake
    dydt[idx["A_brain_ISF"]] += (
        - fluxes["ISF_to_cell"]
    )

    # A_brain_cell:
    # gain from ISF, lose to nucleus and cell-processing loss
    dydt[idx["A_brain_cell"]] += (
        fluxes["ISF_to_cell"]
        - fluxes["cell_to_nuc"]
        - fluxes["cell_loss"]
    )

    # A_brain_nuc:
    # gain from cell, lose by degradation
    dydt[idx["A_brain_nuc"]] += (
        fluxes["cell_to_nuc"]
        - fluxes["deg_v"]
    )

    # mRNA_brain:
    # gain by transcription, lose by degradation
    dydt[idx["mRNA_brain"]] += (
        fluxes["tx"]
        - fluxes["deg_m"]
    )

    # P_brain:
    # gain by translation, lose by degradation
    dydt[idx["P_brain"]] += (
        fluxes["tl"]
        - fluxes["deg_p"]
    )

    # Optional accumulation into cleared pool
    dydt[idx["A_cleared"]] += (
        fluxes["cell_loss"]
        + fluxes["deg_v"]
        + fluxes["deg_m"]
        + fluxes["deg_p"]
    )

    return dydt


def extract_intracellular_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct module-4 intracellular flux trajectories from simulation result.
    """
    flux_history = {
        "ISF_to_cell": [],
        "cell_to_nuc": [],
        "cell_loss": [],
        "deg_v": [],
        "tx": [],
        "deg_m": [],
        "tl": [],
        "deg_p": [],
    }

    for k in range(y_array.shape[1]):
        A_brain_ISF = y_array[idx["A_brain_ISF"], k]
        A_brain_cell = y_array[idx["A_brain_cell"], k]
        A_brain_nuc = y_array[idx["A_brain_nuc"], k]
        mRNA_brain = y_array[idx["mRNA_brain"], k]
        P_brain = y_array[idx["P_brain"], k]

        fluxes = compute_intracellular_fluxes(
            A_brain_ISF=A_brain_ISF,
            A_brain_cell=A_brain_cell,
            A_brain_nuc=A_brain_nuc,
            mRNA_brain=mRNA_brain,
            P_brain=P_brain,
            params=params,
        )

        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history