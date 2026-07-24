import numpy as np
from scipy.integrate import solve_ivp

from models.full_model.state_vector import STATE_ORDER, build_index_map
from models.full_model.rhs_aav import rhs_aav


def build_initial_state(dose: float, config) -> np.ndarray:
    """
    Build the global initial state vector.
    """
    idx = build_index_map(STATE_ORDER)
    y0 = np.zeros(len(STATE_ORDER), dtype=float)

    # Module 1 initial input
    y0[idx["A_dep"]] = dose

    # Module 5 initial substrate pools
    editing_model = config["editing"].get("model", "competitive")
    if editing_model == "apoe_multisite":
        editing = config["editing"]
        y0[idx["S_APOE4"]] = editing["S_APOE4_init"]
        y0[idx["S_puf_off"]] = editing["S_puf_off_init"]
        y0[idx["S_deaminase_bg"]] = editing["S_deaminase_bg_init"]
    else:
        y0[idx["S_on"]] = config["editing"]["S_on_init"]
        y0[idx["S_off"]] = config["editing"]["S_off_init"]

    return y0


def run_simulation(config, dose=1.0, t_end=72.0, dt=0.1):
    """
    Run the current AAV dynamics simulation.
    """
    idx = build_index_map(STATE_ORDER)
    y0 = build_initial_state(dose, config)

    t_eval = np.arange(0.0, t_end + dt, dt)

    sol = solve_ivp(
        fun=lambda t, y: rhs_aav(t, y, config, idx),
        t_span=(0.0, t_end),
        y0=y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-6,
        atol=1e-9,
    )

    return {
        "t": sol.t,
        "y": sol.y,
        "idx": idx,
        "state_order": STATE_ORDER,
        "dose": dose,
        "success": sol.success,
        "message": sol.message,
        "config": config,
    }
