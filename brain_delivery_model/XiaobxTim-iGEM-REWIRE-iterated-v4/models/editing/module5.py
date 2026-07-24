from __future__ import annotations

from typing import Dict
import numpy as np

from models.editing.competitive_editing import (
    competitive_editing_rhs,
    extract_editing_flux_trajectories as extract_competitive_flux_trajectories,
    compute_editing_metrics as compute_competitive_metrics,
)
from models.editing.puf_apobec_mechanistic import (
    puf_apobec_mechanistic_rhs,
    extract_mechanistic_flux_trajectories,
    compute_mechanistic_editing_metrics,
)
from models.editing.apoe_multisite_editing import (
    apoe_multisite_rhs,
    extract_apoe_multisite_flux_trajectories,
    compute_apoe_multisite_metrics,
)


def get_editing_model(params: Dict) -> str:
    """
    Return the selected module-5 model.

    Supported values:
    - competitive: original simplified competitive editing model
    - mechanistic_puf_apobec: explicit binding/unbinding/catalysis model
    - apoe_multisite: APOE112/APOE158-aware v4 model
    """
    return str(params.get("model", "competitive"))


def editing_rhs(t: float, y: np.ndarray, params: Dict, idx: Dict[str, int]) -> np.ndarray:
    model = get_editing_model(params)
    if model == "competitive":
        return competitive_editing_rhs(t, y, params, idx)
    if model == "mechanistic_puf_apobec":
        return puf_apobec_mechanistic_rhs(t, y, params, idx)
    if model == "apoe_multisite":
        return apoe_multisite_rhs(t, y, params, idx)
    raise ValueError(f"Unknown editing model '{model}'.")


def extract_editing_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict,
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    model = get_editing_model(params)
    if model == "competitive":
        return extract_competitive_flux_trajectories(t_array, y_array, params, idx)
    if model == "mechanistic_puf_apobec":
        return extract_mechanistic_flux_trajectories(t_array, y_array, params, idx)
    if model == "apoe_multisite":
        return extract_apoe_multisite_flux_trajectories(t_array, y_array, params, idx)
    raise ValueError(f"Unknown editing model '{model}'.")


def compute_editing_metrics(y: np.ndarray, idx: Dict[str, int], config: Dict) -> Dict[str, np.ndarray]:
    model = get_editing_model(config["editing"])
    if model == "competitive":
        return compute_competitive_metrics(y, idx, config)
    if model == "mechanistic_puf_apobec":
        return compute_mechanistic_editing_metrics(y, idx, config)
    if model == "apoe_multisite":
        return compute_apoe_multisite_metrics(y, idx, config)
    raise ValueError(f"Unknown editing model '{model}'.")
