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
    "E_on",
    "E_off",
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
