from __future__ import annotations

from typing import Dict
import numpy as np
from scipy.optimize import differential_evolution, minimize

from fitting.objective_module1 import objective_module1
from fitting.fit_utils import (
    vector_to_param_dict,
    param_dict_to_vector,
    bounds_dict_to_list,
)

# -----------------------------
# Hard bounds: optimization search space
# -----------------------------
DEFAULT_BOUNDS = {
    "k_lymph": (1e-4, 1.0),
    "k_blood": (1e-4, 1.0),
    "k_deg_loc": (1e-5, 0.1),
    "k_lymph_to_blood": (1e-4, 1.0),
}

# -----------------------------
# Soft range representative initial guess
# -----------------------------
DEFAULT_INIT_PARAMS = {
    "k_lymph": 5e-2,
    "k_blood": 5e-2,
    "k_deg_loc": 1e-3,
    "k_lymph_to_blood": 8e-2,
}


def fit_module1_global(
    dose: float,
    observations: Dict,
    bounds_dict: Dict[str, tuple] | None = None,
    loss_type: str = "log_sse",
    seed: int = 42,
    use_soft_penalty: bool = False,
    penalty_weight: float = 1.0,
) -> Dict:
    """
    Global optimization using differential evolution.
    """
    bounds_dict = DEFAULT_BOUNDS if bounds_dict is None else bounds_dict
    bounds = bounds_dict_to_list(bounds_dict)

    result = differential_evolution(
        func=lambda x: objective_module1(
            x=x,
            dose=dose,
            observations=observations,
            loss_type=loss_type,
            use_soft_penalty=use_soft_penalty,
            penalty_weight=penalty_weight,
        ),
        bounds=bounds,
        seed=seed,
        polish=False,
    )

    return {
        "x": result.x,
        "params": vector_to_param_dict(result.x),
        "loss": float(result.fun),
        "success": bool(result.success),
        "message": result.message,
        "method": "differential_evolution",
    }


def fit_module1_local(
    x0: np.ndarray | None,
    dose: float,
    observations: Dict,
    bounds_dict: Dict[str, tuple] | None = None,
    loss_type: str = "log_sse",
    use_soft_penalty: bool = False,
    penalty_weight: float = 1.0,
) -> Dict:
    """
    Local refinement using L-BFGS-B.
    """
    bounds_dict = DEFAULT_BOUNDS if bounds_dict is None else bounds_dict
    bounds = bounds_dict_to_list(bounds_dict)

    if x0 is None:
        x0 = param_dict_to_vector(DEFAULT_INIT_PARAMS)

    result = minimize(
        fun=lambda x: objective_module1(
            x=x,
            dose=dose,
            observations=observations,
            loss_type=loss_type,
            use_soft_penalty=use_soft_penalty,
            penalty_weight=penalty_weight,
        ),
        x0=x0,
        method="L-BFGS-B",
        bounds=bounds,
    )

    return {
        "x": result.x,
        "params": vector_to_param_dict(result.x),
        "loss": float(result.fun),
        "success": bool(result.success),
        "message": result.message,
        "method": "L-BFGS-B",
    }


def fit_module1_two_stage(
    dose: float,
    observations: Dict,
    bounds_dict: Dict[str, tuple] | None = None,
    loss_type: str = "log_sse",
    seed: int = 42,
    use_soft_penalty: bool = False,
    penalty_weight: float = 1.0,
) -> Dict:
    """
    Two-stage fitting:
    1) global search
    2) local refinement
    """
    global_result = fit_module1_global(
        dose=dose,
        observations=observations,
        bounds_dict=bounds_dict,
        loss_type=loss_type,
        seed=seed,
        use_soft_penalty=use_soft_penalty,
        penalty_weight=penalty_weight,
    )

    local_result = fit_module1_local(
        x0=global_result["x"],
        dose=dose,
        observations=observations,
        bounds_dict=bounds_dict,
        loss_type=loss_type,
        use_soft_penalty=use_soft_penalty,
        penalty_weight=penalty_weight,
    )

    return {
        "global_result": global_result,
        "local_result": local_result,
        "best_params": local_result["params"],
        "best_loss": local_result["loss"],
    }