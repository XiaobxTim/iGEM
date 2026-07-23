from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Any
import pandas as pd
import matplotlib.pyplot as plt


def _ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_dose_scan_table(scan_results: List[Dict[str, Any]], output_dir: str | Path) -> Path:
    out = _ensure_dir(output_dir)
    df = pd.DataFrame(scan_results)
    path = out / "module6_dose_scan.csv"
    df.to_csv(path, index=False)
    return path


def plot_dose_response(scan_results: List[Dict[str, Any]], output_dir: str | Path) -> Dict[str, Path]:
    out = _ensure_dir(output_dir)
    df = pd.DataFrame(scan_results)

    figure_paths = {}

    plot_specs = [
        ("on_target_editing_rate_final", "Dose vs On-target Editing Rate", "dose_vs_on_target.png"),
        ("off_target_burden_final", "Dose vs Off-target Burden", "dose_vs_off_target.png"),
        ("specificity_index_final", "Dose vs Specificity Index", "dose_vs_specificity.png"),
        ("AUC_liver", "Dose vs Liver AUC", "dose_vs_auc_liver.png"),
        ("Cmax_blood", "Dose vs Blood Cmax", "dose_vs_cmax_blood.png"),
        ("P_brain_peak", "Dose vs Brain Expression Peak", "dose_vs_pbrain_peak.png"),
    ]

    for metric_name, title, filename in plot_specs:
        plt.figure(figsize=(7, 4.5))
        plt.plot(df["dose"], df[metric_name], marker="o")
        plt.xlabel("Dose")
        plt.ylabel(metric_name)
        plt.title(title)
        plt.tight_layout()
        path = out / filename
        plt.savefig(path, dpi=160)
        plt.close()
        figure_paths[metric_name] = path

    # Feasibility map
    plt.figure(figsize=(7, 4.5))
    colors = ["green" if x else "red" for x in df["feasible"]]
    plt.scatter(df["dose"], df["on_target_editing_rate_final"], c=colors)
    plt.xlabel("Dose")
    plt.ylabel("on_target_editing_rate_final")
    plt.title("Feasibility Map (green = feasible)")
    plt.tight_layout()
    path = out / "dose_feasibility_map.png"
    plt.savefig(path, dpi=160)
    plt.close()
    figure_paths["feasibility_map"] = path

    return figure_paths