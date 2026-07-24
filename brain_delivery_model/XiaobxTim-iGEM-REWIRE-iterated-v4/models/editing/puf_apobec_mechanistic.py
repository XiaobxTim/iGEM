from __future__ import annotations

from typing import Dict
import numpy as np


REQUIRED_MECHANISTIC_KEYS = (
    "k_on_on",
    "k_off_on",
    "k_cat_on",
    "k_on_off",
    "k_off_off",
    "k_cat_off",
    "k_prod_on",
    "k_deg_on",
    "k_prod_off",
    "k_deg_off",
    "k_loss_on",
    "k_loss_off",
)


def validate_mechanistic_params(params: Dict[str, float]) -> None:
    """
    Validate the explicit PUF-APOBEC editing parameters.

    This module combines two ideas:
    - The GitHub model supplies P_brain, the active editor made after AAV
      delivery, BBB transport, cell uptake, transcription, and translation.
    - The local notebook describes enzyme kinetics: E + S <-> ES -> E + P.

    Here P_brain is treated as the free active editor E. It can bind target
    APOE4 RNA or off-target RNA, form enzyme-substrate complexes, catalyze
    editing, and then be released for reuse.
    """
    for key in REQUIRED_MECHANISTIC_KEYS:
        if key not in params:
            raise KeyError(f"Missing mechanistic editing parameter: {key}")
        if params[key] < 0:
            raise ValueError(f"Mechanistic editing parameter '{key}' must be non-negative.")


def compute_mechanistic_editing_fluxes(
    P_brain: float,
    S_on: float,
    S_off: float,
    ES_on: float,
    ES_off: float,
    E_on: float,
    E_off: float,
    params: Dict[str, float],
) -> Dict[str, float]:
    """
    Compute all reaction fluxes for the combined PUF-APOBEC module.

    Naming convention
    -----------------
    P_brain:
        Free active PUF-APOBEC editor in the brain cell compartment. This is
        produced by module 4 and plays the role of E in the notebook.
    S_on / S_off:
        Unedited on-target APOE4 RNA and accessible off-target RNA pools.
    ES_on / ES_off:
        Explicit enzyme-substrate complexes for the two RNA pools.
    E_on / E_off:
        Edited products accumulated from on-target and off-target catalysis.

    The two substrates compete because both binding fluxes consume the same
    free editor pool P_brain.
    """
    validate_mechanistic_params(params)

    # Numerical guard: LSODA can briefly propose tiny negative values during
    # stiff transitions. Clipping only for flux calculation avoids artificial
    # negative reaction rates while leaving the solver state untouched.
    E_free = max(float(P_brain), 0.0)
    S_on_free = max(float(S_on), 0.0)
    S_off_free = max(float(S_off), 0.0)
    C_on = max(float(ES_on), 0.0)
    C_off = max(float(ES_off), 0.0)

    bind_on = params["k_on_on"] * E_free * S_on_free
    unbind_on = params["k_off_on"] * C_on
    edit_on = params["k_cat_on"] * C_on

    bind_off = params["k_on_off"] * E_free * S_off_free
    unbind_off = params["k_off_off"] * C_off
    edit_off = params["k_cat_off"] * C_off

    return {
        "prod_on": params["k_prod_on"],
        "deg_on": params["k_deg_on"] * S_on_free,
        "bind_on": bind_on,
        "unbind_on": unbind_on,
        "edit_on": edit_on,
        "prod_off": params["k_prod_off"],
        "deg_off": params["k_deg_off"] * S_off_free,
        "bind_off": bind_off,
        "unbind_off": unbind_off,
        "edit_off": edit_off,
        "loss_on": params["k_loss_on"] * max(float(E_on), 0.0),
        "loss_off": params["k_loss_off"] * max(float(E_off), 0.0),
    }


def puf_apobec_mechanistic_rhs(
    t: float,
    y: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> np.ndarray:
    """
    RHS for the mechanistic PUF-APOBEC editing module.

    Compared with the original competitive Michaelis-Menten style module, this
    explicitly tracks transient complexes:

        P_brain + S_on  <->  ES_on   ->  P_brain + E_on
        P_brain + S_off <->  ES_off  ->  P_brain + E_off

    P_brain is also changed by module 4 in the same global RHS. The updates
    below only add the binding/catalysis contribution to the editor pool.
    """
    _ = t
    validate_mechanistic_params(params)

    dydt = np.zeros_like(y)

    P_brain = y[idx["P_brain"]]
    S_on = y[idx["S_on"]]
    S_off = y[idx["S_off"]]
    ES_on = y[idx["ES_on"]]
    ES_off = y[idx["ES_off"]]
    E_on = y[idx["E_on"]]
    E_off = y[idx["E_off"]]

    fluxes = compute_mechanistic_editing_fluxes(
        P_brain=P_brain,
        S_on=S_on,
        S_off=S_off,
        ES_on=ES_on,
        ES_off=ES_off,
        E_on=E_on,
        E_off=E_off,
        params=params,
    )

    # Free active editor. Binding temporarily removes editor from the free
    # pool; unbinding and catalysis release it back, matching the notebook's
    # enzyme recycling term.
    dydt[idx["P_brain"]] += (
        - fluxes["bind_on"]
        - fluxes["bind_off"]
        + fluxes["unbind_on"]
        + fluxes["unbind_off"]
        + fluxes["edit_on"]
        + fluxes["edit_off"]
    )

    # On-target APOE4 RNA: renewed by transcription/turnover, naturally
    # degraded, consumed by binding, and restored if the complex dissociates.
    dydt[idx["S_on"]] += (
        fluxes["prod_on"]
        - fluxes["deg_on"]
        - fluxes["bind_on"]
        + fluxes["unbind_on"]
    )

    # Off-target accessible RNA pool follows the same binding logic.
    dydt[idx["S_off"]] += (
        fluxes["prod_off"]
        - fluxes["deg_off"]
        - fluxes["bind_off"]
        + fluxes["unbind_off"]
    )

    # Enzyme-substrate complexes are formed by binding, then lost by either
    # dissociation or catalytic editing.
    dydt[idx["ES_on"]] += (
        fluxes["bind_on"]
        - fluxes["unbind_on"]
        - fluxes["edit_on"]
    )
    dydt[idx["ES_off"]] += (
        fluxes["bind_off"]
        - fluxes["unbind_off"]
        - fluxes["edit_off"]
    )

    # Edited product pools. E_on is the desired APOE2-like edited RNA output;
    # E_off is the safety burden from off-target edits.
    dydt[idx["E_on"]] += fluxes["edit_on"] - fluxes["loss_on"]
    dydt[idx["E_off"]] += fluxes["edit_off"] - fluxes["loss_off"]

    return dydt


def extract_mechanistic_flux_trajectories(
    t_array: np.ndarray,
    y_array: np.ndarray,
    params: Dict[str, float],
    idx: Dict[str, int],
) -> Dict[str, np.ndarray]:
    """
    Reconstruct mechanistic editing fluxes for plotting and diagnostics.
    """
    _ = t_array
    flux_history = {
        "prod_on": [],
        "deg_on": [],
        "bind_on": [],
        "unbind_on": [],
        "edit_on": [],
        "prod_off": [],
        "deg_off": [],
        "bind_off": [],
        "unbind_off": [],
        "edit_off": [],
        "loss_on": [],
        "loss_off": [],
    }

    for k in range(y_array.shape[1]):
        fluxes = compute_mechanistic_editing_fluxes(
            P_brain=y_array[idx["P_brain"], k],
            S_on=y_array[idx["S_on"], k],
            S_off=y_array[idx["S_off"], k],
            ES_on=y_array[idx["ES_on"], k],
            ES_off=y_array[idx["ES_off"], k],
            E_on=y_array[idx["E_on"], k],
            E_off=y_array[idx["E_off"], k],
            params=params,
        )
        for name in flux_history:
            flux_history[name].append(fluxes[name])

    for name in flux_history:
        flux_history[name] = np.asarray(flux_history[name], dtype=float)

    return flux_history


def compute_mechanistic_editing_metrics(
    y: np.ndarray,
    idx: Dict[str, int],
    config: Dict,
) -> Dict[str, np.ndarray]:
    """
    Compute observable editing metrics.

    The main editing fraction uses edited / (unedited + edited), which matches
    the local notebook. The initial-pool normalized metrics are kept for module
    6 dose optimization compatibility.
    """
    editing_cfg = config["editing"]
    eps = float(editing_cfg.get("specificity_eps", 1e-12))

    S_on = y[idx["S_on"], :]
    S_off = y[idx["S_off"], :]
    E_on = y[idx["E_on"], :]
    E_off = y[idx["E_off"], :]

    S_on_init = float(editing_cfg.get("S_on_init", max(float(S_on[0]), eps)))
    S_off_init = float(editing_cfg.get("S_off_init", max(float(S_off[0]), eps)))

    on_target_editing_fraction = E_on / (S_on + E_on + eps)
    off_target_editing_fraction = E_off / (S_off + E_off + eps)

    return {
        "on_target_editing_rate": E_on / max(S_on_init, eps),
        "off_target_burden": E_off / max(S_off_init, eps),
        "specificity_index": E_on / (E_off + eps),
        "on_target_editing_fraction": on_target_editing_fraction,
        "off_target_editing_fraction": off_target_editing_fraction,
    }
