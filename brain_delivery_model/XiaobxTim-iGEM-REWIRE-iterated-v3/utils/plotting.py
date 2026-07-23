from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import matplotlib.pyplot as plt
import numpy as np

from models.pbpk.lymphatic_absorption import extract_absorption_flux_trajectories
from models.pbpk.organ_distribution import extract_distribution_flux_trajectories
from models.bbb.bbb_transport import extract_bbb_flux_trajectories
from models.intracellular.aav_intracellular import extract_intracellular_flux_trajectories
from models.editing.competitive_editing import (
    extract_editing_flux_trajectories,
    compute_editing_metrics,
)


def _ensure_output_dir(output_dir: str | Path) -> Path:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def _has_state(idx: Dict[str, int], name: str) -> bool:
    return name in idx


def plot_module1_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    plt.figure(figsize=(8, 5))

    for state_name in ["A_dep", "A_lymph", "A_blood"]:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount")
    plt.title("Module 1: Absorption States")
    plt.legend()
    plt.tight_layout()

    path = out / "module1_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path


def plot_module1_fluxes(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    config = results["config"]
    route = config.get("route", "footpad")
    params = config["absorption"][route]

    # 这里用位置参数最稳，避免关键字参数名不一致
    fluxes = extract_absorption_flux_trajectories(t, y, params, idx)

    plt.figure(figsize=(8, 5))
    for name, arr in fluxes.items():
        plt.plot(t, arr, label=name)

    plt.xlabel("Time (h)")
    plt.ylabel("Flux")
    plt.title("Module 1: Absorption Fluxes")
    plt.legend()
    plt.tight_layout()

    path = out / "module1_fluxes.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path


def plot_module2_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    plt.figure(figsize=(9, 5))

    for state_name in ["A_blood", "A_liver", "A_peripheral", "A_brain_blood", "A_cleared"]:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount")
    plt.title("Module 2: Systemic Distribution States")
    plt.legend()
    plt.tight_layout()

    path = out / "module2_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path


def plot_module2_fluxes(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    config = results["config"]
    params = config["distribution"]

    fluxes = extract_distribution_flux_trajectories(
        t,
        y,
        params,
        idx,
    )

    plt.figure(figsize=(10, 6))
    for name, arr in fluxes.items():
        plt.plot(t, arr, label=name)

    plt.xlabel("Time (h)")
    plt.ylabel("Flux")
    plt.title("Module 2: Distribution Fluxes")
    plt.legend()
    plt.tight_layout()

    path = out / "module2_fluxes.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module3_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    plt.figure(figsize=(9, 5))

    for state_name in ["A_brain_blood", "A_brain_EC", "A_brain_endo", "A_brain_ISF"]:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount")
    plt.title("Module 3: BBB Transport States")
    plt.legend()
    plt.tight_layout()

    path = out / "module3_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module3_fluxes(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    config = results["config"]
    params = config["bbb"]

    fluxes = extract_bbb_flux_trajectories(
        t,
        y,
        params,
        idx,
    )

    plt.figure(figsize=(10, 6))
    for name, arr in fluxes.items():
        plt.plot(t, arr, label=name)

    plt.xlabel("Time (h)")
    plt.ylabel("Flux")
    plt.title("Module 3: BBB Transport Fluxes")
    plt.legend()
    plt.tight_layout()

    path = out / "module3_fluxes.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module4_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    plt.figure(figsize=(9, 5))

    for state_name in ["A_brain_ISF", "A_brain_cell", "A_brain_nuc", "mRNA_brain", "P_brain"]:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount / Level")
    plt.title("Module 4: Brain Uptake and Expression States")
    plt.legend()
    plt.tight_layout()

    path = out / "module4_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module4_fluxes(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    config = results["config"]
    params = config["intracellular"]

    fluxes = extract_intracellular_flux_trajectories(
        t,
        y,
        params,
        idx,
    )

    plt.figure(figsize=(10, 6))
    for name, arr in fluxes.items():
        plt.plot(t, arr, label=name)

    plt.xlabel("Time (h)")
    plt.ylabel("Flux")
    plt.title("Module 4: Brain Uptake and Expression Fluxes")
    plt.legend()
    plt.tight_layout()

    path = out / "module4_fluxes.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module5_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    plt.figure(figsize=(9, 5))

    for state_name in ["S_on", "S_off", "E_on", "E_off"]:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount / Level")
    plt.title("Module 5: Competitive Editing States")
    plt.legend()
    plt.tight_layout()

    path = out / "module5_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module5_fluxes(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    config = results["config"]
    params = config["editing"]

    fluxes = extract_editing_flux_trajectories(
        t,
        y,
        params,
        idx,
    )

    plt.figure(figsize=(9, 5))
    for name, arr in fluxes.items():
        plt.plot(t, arr, label=name)

    plt.xlabel("Time (h)")
    plt.ylabel("Flux")
    plt.title("Module 5: Competitive Editing Fluxes")
    plt.legend()
    plt.tight_layout()

    path = out / "module5_fluxes.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_module5_metrics(results: Dict[str, Any], output_dir: str | Path) -> Path:
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]
    config = results["config"]

    metrics = compute_editing_metrics(y=y, idx=idx, config=config)

    plt.figure(figsize=(9, 5))
    plt.plot(t, metrics["on_target_editing_rate"], label="on_target_editing_rate")
    plt.plot(t, metrics["off_target_burden"], label="off_target_burden")
    plt.plot(t, metrics["specificity_index"], label="specificity_index")

    plt.xlabel("Time (h)")
    plt.ylabel("Metric")
    plt.title("Module 5: Editing Metrics")
    plt.legend()
    plt.tight_layout()

    path = out / "module5_metrics.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path

def plot_combined_states(results: Dict[str, Any], output_dir: str | Path) -> Path:
    """
    一张总览图，把模块一和模块二关键状态放一起。
    """
    out = _ensure_output_dir(output_dir)
    t = results["t"]
    y = results["y"]
    idx = results["idx"]

    selected_states = [
        "A_dep",
        "A_lymph",
        "A_blood",
        "A_liver",
        "A_peripheral",
        "A_brain_blood",
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
        "A_cleared",
    ]

    plt.figure(figsize=(10, 6))
    for state_name in selected_states:
        if _has_state(idx, state_name):
            plt.plot(t, y[idx[state_name]], label=state_name)

    plt.xlabel("Time (h)")
    plt.ylabel("Amount")
    plt.title("Key State Trajectories")
    plt.legend()
    plt.tight_layout()

    path = out / "combined_states.png"
    plt.savefig(path, dpi=160)
    plt.close()
    return path


def plot_all(results: Dict[str, Any], output_dir: str | Path) -> Dict[str, Path]:
    figure_paths = {
        "module1_states": plot_module1_states(results, output_dir),
        "module1_fluxes": plot_module1_fluxes(results, output_dir),
        "module2_states": plot_module2_states(results, output_dir),
        "module2_fluxes": plot_module2_fluxes(results, output_dir),
        "module3_states": plot_module3_states(results, output_dir),
        "module3_fluxes": plot_module3_fluxes(results, output_dir),
        "module4_states": plot_module4_states(results, output_dir),
        "module4_fluxes": plot_module4_fluxes(results, output_dir),
        "module5_states": plot_module5_states(results, output_dir),
        "module5_fluxes": plot_module5_fluxes(results, output_dir),
        "module5_metrics": plot_module5_metrics(results, output_dir),
        "combined_states": plot_combined_states(results, output_dir),
    }
    return figure_paths