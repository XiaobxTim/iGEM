from __future__ import annotations

from pathlib import Path
import csv
import argparse
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.config_loader import load_base_config
from utils.io import ensure_dir, save_results
from utils.plotting import plot_all
from models.full_model.simulator import run_simulation
from optimization.metrics import extract_summary_metrics
from optimization.iterated.uncertainty import (
    load_parameter_priors,
    run_uncertainty_scan,
    save_rows,
)
from models.calibration.wetlab_bridge import (
    load_wetlab_observations,
    residual_table,
    weighted_rmse,
)
from models.calibration.parameter_fit import (
    random_search_fit,
    save_fit_config,
    save_fit_results,
)
from models.editing.offtarget_panel import apply_offtarget_panel_to_config
from models.editing.sequence_to_kinetics import (
    apply_design_to_config,
    load_design_table,
    pareto_front,
)


def write_summary(summary: dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])


def write_residuals(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_key_value_table(data: dict, path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["key", "value"])
        for key, value in data.items():
            writer.writerow([key, value])


def save_dict_rows(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def utility(summary: dict, weights: dict | None = None) -> float:
    weights = weights or {
        "apoe3": 1.0,
        "apoe2_risk": 0.7,
        "off": 1.2,
        "local": 0.6,
        "puf": 1.0,
        "deam": 1.0,
        "burden": 0.15,
        "late_time": 0.001,
        "reference_time_h": 48.0,
    }
    sample_time = float(summary.get("_sample_time_h", weights["reference_time_h"]))
    late_time_penalty = weights["late_time"] * max(sample_time - weights["reference_time_h"], 0.0)
    return (
        weights["apoe3"] * float(summary.get("apoe3_like_fraction_final", 0.0))
        - weights["apoe2_risk"] * float(summary.get("ldlr_binding_risk_proxy_final", 0.0))
        - weights["off"] * float(summary.get("off_target_burden_final", 0.0))
        - weights["local"] * float(summary.get("local_bystander_burden_final", 0.0))
        - weights["puf"] * float(summary.get("puf_mismatch_burden_final", 0.0))
        - weights["deam"] * float(summary.get("deaminase_background_burden_final", 0.0))
        - weights["burden"] * float(summary.get("P_brain_peak", 0.0))
        - late_time_penalty
    )


def run_design_screen(base_config: dict, design_path: Path, dose: float, dt: float) -> list[dict]:
    dose_levels = {
        "low": 0.5 * dose,
        "medium": dose,
        "high": 2.0 * dose,
    }
    sample_times = [24.0, 48.0, 72.0]
    rows = []
    for design in load_design_table(design_path):
        design_config, modifiers = apply_design_to_config(base_config, design)
        for dose_label, design_dose in dose_levels.items():
            for sample_time in sample_times:
                sim = run_simulation(
                    config=design_config,
                    dose=design_dose,
                    t_end=sample_time,
                    dt=dt,
                )
                summary = extract_summary_metrics(sim)
                summary["_sample_time_h"] = sample_time
                rows.append(
                    {
                        "design_id": design["design_id"],
                        "editor_type": design["editor_type"],
                        "dose_level": dose_label,
                        "dose": design_dose,
                        "sample_time_h": sample_time,
                        "utility": utility(summary),
                        **modifiers,
                        **summary,
                    }
                )
    rows.sort(key=lambda row: row["utility"], reverse=True)
    return rows


def write_design_recommendation(rows: list[dict], path: Path) -> None:
    front = pareto_front(rows)
    front_ids = {(
        row["design_id"],
        row["dose_level"],
        row["sample_time_h"],
    ) for row in front}
    with open(path, "w", encoding="utf-8") as f:
        f.write("# v4 design recommendation\n\n")
        if rows:
            best = rows[0]
            f.write("## Recommended design\n\n")
            f.write(f"- Design: `{best['design_id']}`\n")
            f.write(f"- Editor: `{best['editor_type']}`\n")
            f.write(f"- Plasmid/expression level: `{best['dose_level']}`\n")
            f.write(f"- Sampling time: `{best['sample_time_h']} h`\n")
            f.write(f"- Utility: `{best['utility']:.4g}`\n")
            f.write(f"- Expected APOE4 fraction: `{best.get('apoe4_fraction_final', 0.0):.4g}`\n")
            f.write(f"- Expected APOE3-like fraction: `{best.get('apoe3_like_fraction_final', 0.0):.4g}`\n")
            f.write(f"- Expected APOE2-like risk proxy: `{best.get('apoe2_like_fraction_final', 0.0):.4g}`\n")
            f.write(f"- Off-target burden: `{best.get('off_target_burden_final', 0.0):.4g}`\n\n")
        f.write("## Top design rows\n\n")
        f.write("| rank | design | editor | level | time_h | utility | APOE3-like | APOE2-risk | off-target | Pareto |\n")
        f.write("|---:|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for rank, row in enumerate(rows[:12], start=1):
            key = (row["design_id"], row["dose_level"], row["sample_time_h"])
            is_front = "yes" if key in front_ids else "no"
            f.write(
                f"| {rank} | `{row['design_id']}` | `{row['editor_type']}` | `{row['dose_level']}` | "
                f"{row['sample_time_h']:.0f} | {row['utility']:.4g} | "
                f"{row.get('apoe3_like_fraction_final', 0.0):.4g} | "
                f"{row.get('apoe2_like_fraction_final', 0.0):.4g} | "
                f"{row.get('off_target_burden_final', 0.0):.4g} | {is_front} |\n"
            )


def _pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 3:
        return 0.0
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    x_var = sum((x - x_mean) ** 2 for x in xs)
    y_var = sum((y - y_mean) ** 2 for y in ys)
    if x_var <= 0 or y_var <= 0:
        return 0.0
    return numerator / ((x_var * y_var) ** 0.5)


def write_experimental_design_recommendations(rows: list[dict], path: Path) -> None:
    """
    Rank uncertain parameters by influence on efficacy/safety metrics.

    This is deliberately simple and transparent for iGEM documentation. It is
    not a replacement for global sensitivity analysis, but it is enough to make
    the next wet-lab round model-informed.
    """
    metrics = [
        "on_target_editing_rate_final",
        "off_target_burden_final",
        "specificity_index_final",
        "P_brain_peak",
        "AUC_liver",
        "Cmax_blood",
    ]
    parameter_keys = [
        key
        for key in rows[0].keys()
        if "." in key and key not in metrics
    ]

    rankings = []
    for parameter in parameter_keys:
        xs = [float(row[parameter]) for row in rows if row.get(parameter) not in ("", None)]
        if len(xs) != len(rows):
            continue
        score = 0.0
        influences = {}
        for metric in metrics:
            ys = [float(row[metric]) for row in rows]
            corr = _pearson(xs, ys)
            influences[metric] = corr
            score += abs(corr)
        rankings.append((score, parameter, influences))
    rankings.sort(reverse=True)

    assay_map = {
        "distribution.k_blood_to_liver": "liver vector-genome qPCR/ddPCR time course",
        "distribution.k_blood_to_brain": "brain vector-genome qPCR/ddPCR and tissue fractionation",
        "bbb.k_brainblood_to_EC": "brain endothelial versus parenchymal vector localization",
        "bbb.k_endo_to_ISF": "brain parenchymal delivery assay after perfusion",
        "intracellular.k_tx": "PUF-APOBEC mRNA qPCR time course",
        "intracellular.k_deg_m": "actinomycin-D or washout mRNA stability assay",
        "intracellular.k_tl": "tagged editor fluorescence or Western blot time course",
        "intracellular.k_deg_p": "cycloheximide chase or pulse-chase editor stability assay",
        "editing.k_on_on": "EMSA/SPR or early-time on-target editing kinetics",
        "editing.k_off_on": "EMSA/SPR dissociation or editing washout experiment",
        "editing.k_cat_on": "amplicon sequencing editing time course at fixed editor level",
        "editing.k_on_off": "RNA-seq/amplicon off-target panel binding/editing screen",
        "editing.k_cat_off": "off-target amplicon time course for top predicted sites",
    }

    with open(path, "w", encoding="utf-8") as f:
        f.write("# Model-informed next wet-lab experiments\n\n")
        f.write("Parameters are ranked by summed absolute Pearson correlation across key decision metrics.\n\n")
        f.write("| rank | parameter | score | recommended experiment |\n")
        f.write("|---:|---|---:|---|\n")
        for rank, (score, parameter, _influences) in enumerate(rankings[:10], start=1):
            experiment = assay_map.get(parameter, "targeted calibration experiment")
            f.write(f"| {rank} | `{parameter}` | {score:.3f} | {experiment} |\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the iterated wetlab-integrated model analysis.")
    parser.add_argument("--dose", type=float, default=1.0)
    parser.add_argument("--t_end", type=float, default=72.0)
    parser.add_argument("--dt", type=float, default=0.2)
    parser.add_argument("--uncertainty_samples", type=int, default=24)
    parser.add_argument("--fit_samples", type=int, default=24)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument(
        "--observations",
        type=str,
        default="wetlab/templates/assay_observations_template.csv",
        help="CSV with wet-lab observations or placeholder template values.",
    )
    parser.add_argument(
        "--offtarget_panel",
        type=str,
        default="wetlab/templates/offtarget_candidate_panel_template.csv",
        help="Candidate off-target panel to aggregate into S_off/k_on_off/k_cat_off.",
    )
    parser.add_argument(
        "--design_table",
        type=str,
        default="wetlab/templates/puf_design_candidates_template.csv",
        help="PUF-deaminase design candidate table for v4 Pareto/utility screening.",
    )
    parser.add_argument("--output_dir", type=str, default="reports/iterated_analysis")
    args = parser.parse_args()

    project_root = PROJECT_ROOT
    config = load_base_config(project_root)
    out = ensure_dir(project_root / args.output_dir)

    panel_path = project_root / args.offtarget_panel
    if panel_path.exists():
        config, panel_summary = apply_offtarget_panel_to_config(config, panel_path)
        write_key_value_table(panel_summary, out / "offtarget_panel_summary.csv")

    design_path = project_root / args.design_table
    if design_path.exists():
        design_rows = run_design_screen(config, design_path, args.dose, args.dt)
        save_dict_rows(design_rows, out / "v4_design_screen.csv")
        write_design_recommendation(design_rows, out / "v4_design_recommendation.md")

    baseline = run_simulation(config=config, dose=args.dose, t_end=args.t_end, dt=args.dt)
    save_results(baseline, out / "baseline")
    plot_all(baseline, out / "baseline")

    summary = extract_summary_metrics(baseline)
    write_summary(summary, out / "baseline_summary.csv")

    obs_path = project_root / args.observations
    observations = []
    if obs_path.exists():
        observations = load_wetlab_observations(obs_path)
        rows = residual_table(baseline, observations)
        write_residuals(rows, out / "wetlab_residuals.csv")
        with open(out / "wetlab_objective.txt", "w", encoding="utf-8") as f:
            f.write(f"weighted_rmse={weighted_rmse(baseline, observations):.6g}\n")

    priors = load_parameter_priors(project_root / "config" / "parameter_provenance.yaml")
    if observations and args.fit_samples > 0:
        fitted_config, fit_rows = random_search_fit(
            base_config=config,
            priors=priors,
            observations=observations,
            n_samples=args.fit_samples,
            dose=args.dose,
            t_end=args.t_end,
            dt=args.dt,
            seed=args.seed + 101,
        )
        save_fit_results(fit_rows, out / "parameter_fit_results.csv")
        save_fit_config(fitted_config, out / "best_fit_config.yaml")

        fitted = run_simulation(config=fitted_config, dose=args.dose, t_end=args.t_end, dt=args.dt)
        save_results(fitted, out / "fitted")
        plot_all(fitted, out / "fitted")
        fitted_summary = extract_summary_metrics(fitted)
        write_summary(fitted_summary, out / "fitted_summary.csv")

    uncertainty_rows = list(run_uncertainty_scan(
        base_config=config,
        priors=priors,
        n_samples=args.uncertainty_samples,
        dose=args.dose,
        t_end=args.t_end,
        dt=args.dt,
        seed=args.seed,
    ))
    save_rows(uncertainty_rows, out / "uncertainty_samples.csv")
    if uncertainty_rows:
        write_experimental_design_recommendations(
            uncertainty_rows,
            out / "experimental_design_recommendations.md",
        )

    print("Iterated analysis finished.")
    print(f"Output directory: {out}")
    print(f"Baseline on-target editing final: {summary['on_target_editing_rate_final']:.6g}")
    print(f"Baseline off-target burden final: {summary['off_target_burden_final']:.6g}")
    print(f"Baseline specificity index final: {summary['specificity_index_final']:.6g}")


if __name__ == "__main__":
    main()
