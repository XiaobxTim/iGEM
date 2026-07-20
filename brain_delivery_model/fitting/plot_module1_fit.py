from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import matplotlib.pyplot as plt

from fitting.simulate_module1 import simulate_module1
from models.pbpk.lymphatic_absorption import extract_absorption_flux_trajectories


def _ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_module1_fit_observations(
    best_params: Dict[str, float],
    dose: float,
    observations: Dict,
    output_dir: str | Path,
    t_end: Optional[float] = None,
    n_points: int = 500,
    filename: str = "module1_fit_observations.png",
) -> Path:
    """
    Plot fitted trajectories against observed data points.

    Parameters
    ----------
    best_params : dict
        Best-fit parameter dictionary.
    dose : float
        Initial dose.
    observations : dict
        Observation dictionary.
    output_dir : str or Path
        Directory to save figure.
    t_end : float, optional
        End time for plotting. If None, inferred from observations.
    n_points : int
        Number of simulation points.
    filename : str
        Output figure filename.

    Returns
    -------
    Path
        Saved figure path.
    """
    output_dir = _ensure_dir(output_dir)

    all_times = []
    for obs in observations.values():
        all_times.extend(list(obs["t"]))
    if t_end is None:
        t_end = max(all_times)

    t_eval = np.linspace(0.0, float(t_end), n_points)
    sim = simulate_module1(best_params, dose=dose, t_eval=t_eval)

    idx = sim["idx"]

    observed_states = list(observations.keys())
    n_states = len(observed_states)

    fig, axes = plt.subplots(n_states, 1, figsize=(8, 4 * n_states), squeeze=False)

    for ax, state_name in zip(axes.flat, observed_states):
        y_sim = sim["y"][idx[state_name], :]
        ax.plot(sim["t"], y_sim, label=f"Fit: {state_name}", linewidth=2)

        t_obs = np.asarray(observations[state_name]["t"], dtype=float)
        y_obs = np.asarray(observations[state_name]["y"], dtype=float)
        ax.scatter(t_obs, y_obs, label=f"Obs: {state_name}", zorder=5)

        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Amount / Signal")
        ax.set_title(f"{state_name}: fitted trajectory vs observations")
        ax.grid(True, alpha=0.3)
        ax.legend()

    fig.tight_layout()
    save_path = output_dir / filename
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_module1_fit_states(
    best_params: Dict[str, float],
    dose: float,
    output_dir: str | Path,
    t_end: float = 72.0,
    n_points: int = 500,
    filename: str = "module1_fit_states.png",
) -> Path:
    """
    Plot core module 1 state trajectories.
    """
    output_dir = _ensure_dir(output_dir)
    t_eval = np.linspace(0.0, float(t_end), n_points)
    sim = simulate_module1(best_params, dose=dose, t_eval=t_eval)
    idx = sim["idx"]

    fig, ax = plt.subplots(figsize=(8, 5))
    for state_name in ["A_dep", "A_lymph", "A_blood"]:
        ax.plot(sim["t"], sim["y"][idx[state_name], :], label=state_name, linewidth=2)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Amount")
    ax.set_title("Module 1 fitted state trajectories")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()
    save_path = output_dir / filename
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return save_path


def plot_module1_fit_fluxes(
    best_params: Dict[str, float],
    dose: float,
    output_dir: str | Path,
    t_end: float = 72.0,
    n_points: int = 500,
    filename: str = "module1_fit_fluxes.png",
) -> Path:
    """
    Plot absorption flux trajectories under fitted parameters.
    """
    output_dir = _ensure_dir(output_dir)
    t_eval = np.linspace(0.0, float(t_end), n_points)
    sim = simulate_module1(best_params, dose=dose, t_eval=t_eval)
    idx = sim["idx"]

    fluxes = extract_absorption_flux_trajectories(
        sim["t"],
        sim["y"],
        best_params,
        idx,
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, values in fluxes.items():
        ax.plot(sim["t"], values, label=name, linewidth=2)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Flux")
    ax.set_title("Module 1 fitted absorption fluxes")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()
    save_path = output_dir / filename
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return save_path