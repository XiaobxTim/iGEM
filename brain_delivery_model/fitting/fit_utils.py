from __future__ import annotations

from typing import Dict, List
import numpy as np


MODULE1_PARAM_ORDER = [
    "k_lymph",
    "k_blood",
    "k_deg_loc",
    "k_lymph_to_blood",
]


def vector_to_param_dict(x: np.ndarray) -> Dict[str, float]:
    """
    Convert optimization vector to parameter dictionary.
    """
    return {name: float(val) for name, val in zip(MODULE1_PARAM_ORDER, x)}


def param_dict_to_vector(params: Dict[str, float]) -> np.ndarray:
    """
    Convert parameter dictionary to optimization vector.
    """
    return np.asarray([params[name] for name in MODULE1_PARAM_ORDER], dtype=float)


def bounds_dict_to_list(bounds_dict: Dict[str, tuple]) -> List[tuple]:
    """
    Convert bounds dict to ordered list for optimizers.
    """
    return [bounds_dict[name] for name in MODULE1_PARAM_ORDER]