from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from models.editing.apoe_multisite_editing import apoe_multisite_rhs
from models.full_model.state_vector import IDX, STATE_ORDER


STATE_NAMES = (
    "A_cell",
    "A_nucleus",
    "editor_mRNA",
    "editor_protein",
    "S_APOE4",
    "S_APOE3_like",
    "S_APOE2_like",
    "S_APOE158_only",
    "C_APOE112",
    "C_APOE158",
    "C_APOE158_after112",
    "C_APOE112_after158",
    "B_local_bystander",
    "S_puf_off",
    "C_puf_off",
    "E_puf_off",
    "S_deaminase_bg",
    "E_deaminase_bg",
)

LOCAL_TO_GLOBAL = {
    "A_cell": "A_brain_cell",
    "A_nucleus": "A_brain_nuc",
    "editor_mRNA": "mRNA_brain",
    "editor_protein": "P_brain",
    **{name: name for name in STATE_NAMES[4:]},
}


@dataclass(frozen=True)
class IntracellularTrajectory:
    time_min: np.ndarray
    state_names: tuple[str, ...]
    values: np.ndarray


def _hour_parameters(minute_parameters: dict) -> dict:
    restored: dict[str, object] = {}
    for key, value in minute_parameters.items():
        if key.endswith("_per_min"):
            restored[key.removesuffix("_per_min")] = float(value) * 60.0
        else:
            restored[key] = value
    return restored


def make_initial_state(parameters: dict, apoe_scale: float) -> np.ndarray:
    if not np.isfinite(apoe_scale) or apoe_scale < 0.0:
        raise ValueError("apoe_scale must be finite and non-negative")
    editing = parameters["editing"]
    state = np.zeros(len(STATE_NAMES), dtype=float)
    state[STATE_NAMES.index("S_APOE4")] = float(editing["S_APOE4_init"]) * apoe_scale
    state[STATE_NAMES.index("S_puf_off")] = float(editing["S_puf_off_init"])
    state[STATE_NAMES.index("S_deaminase_bg")] = float(editing["S_deaminase_bg_init"])
    return state


def intracellular_rhs(
    state: np.ndarray,
    parameters: dict,
    uptake_amount_per_min: float,
    apoe_scale: float,
) -> np.ndarray:
    if uptake_amount_per_min < 0.0 or not np.isfinite(uptake_amount_per_min):
        raise ValueError("uptake must be finite and non-negative")
    intracellular = parameters["intracellular"]
    derivative = np.zeros_like(state)
    local = {name: index for index, name in enumerate(STATE_NAMES)}

    a_cell = state[local["A_cell"]]
    a_nucleus = state[local["A_nucleus"]]
    mrna = state[local["editor_mRNA"]]
    protein = state[local["editor_protein"]]
    derivative[local["A_cell"]] = uptake_amount_per_min - (
        intracellular["k_cell_to_nuc_per_min"] + intracellular["k_cell_loss_per_min"]
    ) * a_cell
    derivative[local["A_nucleus"]] = (
        intracellular["k_cell_to_nuc_per_min"] * a_cell
        - intracellular["k_deg_v_per_min"] * a_nucleus
    )
    derivative[local["editor_mRNA"]] = (
        intracellular["k_tx_per_min"] * a_nucleus
        - intracellular["k_deg_m_per_min"] * mrna
    )
    derivative[local["editor_protein"]] = (
        intracellular["k_tl_per_min"] * mrna
        - intracellular["k_deg_p_per_min"] * protein
    )

    global_state = np.zeros(len(STATE_ORDER), dtype=float)
    for local_name, global_name in LOCAL_TO_GLOBAL.items():
        global_state[IDX[global_name]] = state[local[local_name]]
    editing_hour = _hour_parameters(parameters["editing"])
    editing_hour["k_prod_apoe"] = float(editing_hour["k_prod_apoe"]) * apoe_scale
    editing_derivative = apoe_multisite_rhs(0.0, global_state, editing_hour, IDX) / 60.0
    for local_name, global_name in LOCAL_TO_GLOBAL.items():
        if local_name in {"A_cell", "A_nucleus", "editor_mRNA"}:
            continue
        derivative[local[local_name]] += editing_derivative[IDX[global_name]]
    return derivative


def rk4_step(
    state: np.ndarray,
    dt_min: float,
    parameters: dict,
    uptake_amount_per_min: float,
    apoe_scale: float,
) -> np.ndarray:
    rhs = lambda value: intracellular_rhs(
        value, parameters, uptake_amount_per_min, apoe_scale
    )
    k1 = rhs(state)
    k2 = rhs(state + 0.5 * dt_min * k1)
    k3 = rhs(state + 0.5 * dt_min * k2)
    k4 = rhs(state + dt_min * k3)
    updated = state + dt_min * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    if not np.all(np.isfinite(updated)) or np.min(updated) < -1e-9:
        raise FloatingPointError("intracellular integration produced an invalid state")
    return np.maximum(updated, 0.0)


def simulate_intracellular(
    parameters: dict,
    minutes: float,
    dt_min: float,
    uptake_amount_per_min: float,
    apoe_scale: float,
) -> IntracellularTrajectory:
    if minutes <= 0.0 or dt_min <= 0.0 or minutes / dt_min % 1.0 != 0.0:
        raise ValueError("minutes must be a positive integer multiple of dt_min")
    state = make_initial_state(parameters, apoe_scale)
    step_count = int(round(minutes / dt_min))
    values = [state.copy()]
    for _ in range(step_count):
        state = rk4_step(
            state,
            dt_min,
            parameters,
            uptake_amount_per_min,
            apoe_scale,
        )
        values.append(state.copy())
    return IntracellularTrajectory(
        time_min=np.arange(step_count + 1, dtype=float) * dt_min,
        state_names=STATE_NAMES,
        values=np.asarray(values),
    )

