from __future__ import annotations

from typing import Dict, List, Tuple
import numpy as np

STATE_ORDER: List[str] = [
    "A_dep",
    "A_lymph",
    "A_blood",
    "A_liver",
    "A_peripheral",
    "A_brain_blood",
    "A_cleared",
    "A_brain_EC",
    "A_brain_endo",
    "A_brain_ISF",
    "A_brain_cell",
    "A_brain_nuc",
    "mRNA_brain",
    "P_brain",
    "S_on",
    "S_off",
    "ES_on",
    "ES_off",
    "E_on",
    "E_off",
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
]


def build_index_map(state_order: List[str] | None = None) -> Dict[str, int]:
    """Build a name -> index mapping for the global state vector."""
    order = state_order or STATE_ORDER
    return {name: i for i, name in enumerate(order)}


IDX: Dict[str, int] = build_index_map()


def make_initial_state(dose: float) -> np.ndarray:
    """
    Create the default initial condition vector.

    Parameters
    ----------
    dose : float
        Initial administered AAV amount placed in the depot compartment.
    """
    y0 = np.zeros(len(STATE_ORDER), dtype=float)
    y0[IDX["A_dep"]] = dose
    return y0


def as_state_dict(y: np.ndarray, idx: Dict[str, int] | None = None) -> Dict[str, float]:
    """Convert a state vector into a readable dictionary."""
    index_map = idx or IDX
    return {name: float(y[i]) for name, i in index_map.items()}


__all__ = [
    "STATE_ORDER",
    "IDX",
    "build_index_map",
    "make_initial_state",
    "as_state_dict",
]
