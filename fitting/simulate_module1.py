from __future__ import annotations

from typing import Dict, Optional
import numpy as np
from scipy.integrate import solve_ivp

from models.pbpk.lymphatic_absorption import lymphatic_absorption_rhs
from models.full_model.state_vector import STATE_ORDER, build_index_map


def build_module1_initial_state(dose: float) -> np.ndarray:
    """
    Build initial state vector for module 1 only.

    Parameters
    ----------
    dose : float
        Initial dose loaded into depot compartment.

    Returns
    -------
    np.ndarray
        Full state vector with only A_dep initialized.
    """
    idx = build_index_map(STATE_ORDER)
    y0 = np.zeros(len(STATE_ORDER), dtype=float)
    y0[idx["A_dep"]] = dose
    return y0


def simulate_module1(
    params: Dict[str, float],
    dose: float,
    t_eval: np.ndarray,
    method: str = "LSODA",
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> Dict:
    """
    Simulate module 1 (absorption only).

    Parameters
    ----------
    params : dict
        Module 1 parameter dictionary:
        - k_lymph
        - k_blood
        - k_deg_loc
        - k_lymph_to_blood
    dose : float
        Initial dose.
    t_eval : np.ndarray
        Time points for evaluation.
    method : str
        ODE solver method.
    rtol, atol : float
        Solver tolerances.

    Returns
    -------
    dict
        {
            "t": time array,
            "y": state trajectories,
            "idx": state index map,
            "success": bool,
            "message": str
        }
    """
    idx = build_index_map(STATE_ORDER)
    y0 = build_module1_initial_state(dose)

    def rhs(t, y):
        return lymphatic_absorption_rhs(t, y, params, idx)

    sol = solve_ivp(
        rhs,
        (float(t_eval[0]), float(t_eval[-1])),
        y0,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
    )

    return {
        "t": sol.t,
        "y": sol.y,
        "idx": idx,
        "success": sol.success,
        "message": sol.message,
    }