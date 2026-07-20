from __future__ import annotations

from typing import Dict
import numpy as np

from fitting.simulate_module1 import simulate_module1
from fitting.fit_utils import vector_to_param_dict


# -----------------------------
# Soft ranges: preferred parameter region
# Not hard constraints
# -----------------------------
SOFT_RANGES = {
    "k_lymph": (1e-3, 3e-1),
    "k_blood": (1e-3, 3e-1),
    "k_deg_loc": (1e-4, 1e-2),
    "k_lymph_to_blood": (1e-2, 3e-1),
}


def interpolate_state(sim_result: Dict, state_name: str, t_obs: np.ndarray) -> np.ndarray:
    """
    Interpolate simulated state trajectory at observation times.
    """
    idx = sim_result["idx"]
    t_sim = sim_result["t"]
    y_sim = sim_result["y"][idx[state_name], :]
    return np.interp(t_obs, t_sim, y_sim)


def squared_error(y_pred: np.ndarray, y_obs: np.ndarray) -> float:
    return float(np.sum((y_pred - y_obs) ** 2))


def log_squared_error(y_pred: np.ndarray, y_obs: np.ndarray, eps: float = 1e-12) -> float:
    return float(np.sum((np.log(y_pred + eps) - np.log(y_obs + eps)) ** 2))


def soft_range_penalty(
    params: Dict[str, float],
    soft_ranges: Dict[str, tuple] = SOFT_RANGES,
    penalty_weight: float = 1.0,
) -> float:
    """
    Quadratic penalty outside soft preferred ranges.
    No penalty inside range.
    """
    penalty = 0.0
    for name, (low, high) in soft_ranges.items():
        value = params[name]
        if value < low:
            penalty += penalty_weight * (low - value) ** 2
        elif value > high:
            penalty += penalty_weight * (value - high) ** 2
    return float(penalty)


def objective_module1(
    x: np.ndarray,
    dose: float,
    observations: Dict,
    loss_type: str = "log_sse",
    use_soft_penalty: bool = False,
    penalty_weight: float = 1.0,
) -> float:
    """
    Objective function for module 1 fitting.

    Parameters
    ----------
    x : np.ndarray
        Optimization vector in order:
        [k_lymph, k_blood, k_deg_loc, k_lymph_to_blood]
    dose : float
        Initial dose.
    observations : dict
        Observation dictionary, e.g.
        {
            "A_blood": {"t": ..., "y": ..., "weight": 1.0},
            "A_dep": {"t": ..., "y": ..., "weight": 0.5},
        }
    loss_type : str
        "sse" or "log_sse"
    use_soft_penalty : bool
        Whether to add soft-range penalty.
    penalty_weight : float
        Weight of soft-range penalty.

    Returns
    -------
    float
        Total weighted loss.
    """
    params = vector_to_param_dict(x)

    # Build a combined simulation grid that covers all observation times
    all_times = []
    for obs in observations.values():
        all_times.extend(list(obs["t"]))
    t_eval = np.unique(np.asarray(all_times, dtype=float))

    sim_result = simulate_module1(params=params, dose=dose, t_eval=t_eval)

    if not sim_result["success"]:
        return 1e30

    total_loss = 0.0

    for state_name, obs in observations.items():
        t_obs = np.asarray(obs["t"], dtype=float)
        y_obs = np.asarray(obs["y"], dtype=float)
        weight = float(obs.get("weight", 1.0))

        y_pred = interpolate_state(sim_result, state_name, t_obs)

        if loss_type == "sse":
            loss = squared_error(y_pred, y_obs)
        elif loss_type == "log_sse":
            loss = log_squared_error(y_pred, y_obs)
        else:
            raise ValueError(f"Unsupported loss_type: {loss_type}")

        total_loss += weight * loss

    if use_soft_penalty:
        total_loss += soft_range_penalty(
            params=params,
            soft_ranges=SOFT_RANGES,
            penalty_weight=penalty_weight,
        )

    return float(total_loss)