from __future__ import annotations

from pathlib import Path
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from utils.config_loader import load_base_config
from models.full_model.simulator import run_simulation
from optimization.metrics import extract_summary_metrics


PARAMETER_MAP = {
    "k_lymph": ("absorption", "k_lymph"),
    "k_blood": ("absorption", "k_blood"),
    "k_deg_loc": ("absorption", "k_deg_loc"),

    "k_blood_to_liver": ("distribution", "k_blood_to_liver"),
    "k_blood_to_brain": ("distribution", "k_blood_to_brain"),
    "k_clear_blood": ("distribution", "k_clear_blood"),

    "k_endo_to_ISF": ("bbb", "k_endo_to_ISF"),
    "k_endo_loss": ("bbb", "k_endo_loss"),

    "k_ISF_to_cell": ("intracellular", "k_ISF_to_cell"),
    "k_cell_to_nuc": ("intracellular", "k_cell_to_nuc"),
    "k_tx": ("intracellular", "k_tx"),
    "k_tl": ("intracellular", "k_tl"),

    "k_cat_on": ("editing", "k_cat_on"),
    "k_cat_off": ("editing", "k_cat_off"),
}


OUTPUT_METRICS = [
    "P_brain_peak",
    "on_target_editing_rate_final",
    "off_target_burden_final",
    "specificity_index_final",
    "AUC_liver",
    "Cmax_blood",
]


def get_param_block(config: dict, block_name: str) -> dict:
    if block_name == "absorption":
        route = config["route"]
        return config["absorption"][route]
    return config[block_name]


def set_parameter(config: dict, param_name: str, scale: float) -> tuple[float, float]:
    block_name, key = PARAMETER_MAP[param_name]
    block = get_param_block(config, block_name)
    old_value = block[key]
    block[key] = old_value * scale
    return old_value, block[key]


def run_model_and_extract(config: dict, dose: float, t_end: float, dt: float) -> dict:
    results = run_simulation(config=config, dose=dose, t_end=t_end, dt=dt)
    return extract_summary_metrics(results)


def compute_normalized_sensitivity(y_base: float, y_pert: float, p_base: float, p_pert: float) -> float:
    if y_base == 0 or p_base == 0:
        return np.nan
    return ((y_pert - y_base) / y_base) / ((p_pert - p_base) / p_base)


def main():
    # 因为脚本在项目根目录，所以用 .parent
    project_root = Path(__file__).resolve().parent
    config = load_base_config(project_root)

    dose = 1.0
    t_end = 168.0
    dt = 0.2

    output_dir = project_root / "outputs" / "sensitivity_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    base_metrics = run_model_and_extract(copy.deepcopy(config), dose=dose, t_end=t_end, dt=dt)

    rows = []

    for param_name in PARAMETER_MAP:
        for scale in [0.8, 1.2]:
            cfg = copy.deepcopy(config)
            p_base, p_pert = set_parameter(cfg, param_name, scale)
            pert_metrics = run_model_and_extract(cfg, dose=dose, t_end=t_end, dt=dt)

            row = {
                "parameter": param_name,
                "scale": scale,
                "p_base": p_base,
                "p_pert": p_pert,
            }

            for metric in OUTPUT_METRICS:
                row[f"{metric}_base"] = base_metrics[metric]
                row[f"{metric}_pert"] = pert_metrics[metric]
                row[f"{metric}_sens"] = compute_normalized_sensitivity(
                    base_metrics[metric],
                    pert_metrics[metric],
                    p_base,
                    p_pert,
                )

            rows.append(row)

    df = pd.DataFrame(rows)
    csv_path = output_dir / "local_sensitivity_results.csv"
    df.to_csv(csv_path, index=False)

    # 画图：用 +20% 扰动的归一化敏感度
    df_plot = df[df["scale"] == 1.2].copy()

    for metric in OUTPUT_METRICS:
        plot_df = df_plot[["parameter", f"{metric}_sens"]].copy()
        plot_df["abs_sens"] = plot_df[f"{metric}_sens"].abs()
        plot_df = plot_df.sort_values("abs_sens", ascending=False)

        plt.figure(figsize=(9, 5))
        plt.bar(plot_df["parameter"], plot_df[f"{metric}_sens"])
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Normalized sensitivity")
        plt.title(f"Sensitivity of {metric}")
        plt.tight_layout()
        plt.savefig(output_dir / f"sensitivity_{metric}.png", dpi=160)
        plt.close()

    print(f"Sensitivity analysis finished. Results saved to: {output_dir}")
    print(f"CSV: {csv_path}")


if __name__ == "__main__":
    main()