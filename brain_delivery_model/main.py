from __future__ import annotations

import argparse
from pathlib import Path

from utils.config_loader import load_base_config
from utils.io import save_results, ensure_dir
from utils.plotting import plot_all
from models.full_model.simulator import run_simulation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="AAV dynamics model: single simulation or module 6 dose optimization."
    )

    parser.add_argument(
        "--route",
        type=str,
        default=None,
        choices=["footpad", "iv", "im"],
        help="Administration route."
    )
    parser.add_argument(
        "--dose",
        type=float,
        default=1.0,
        help="Initial administered AAV amount for single simulation."
    )
    parser.add_argument(
        "--t_end",
        type=float,
        default=None,
        help="Simulation end time in hours."
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="Simulation output time step."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="outputs/run_full_model",
        help="Directory for single-simulation figures and saved results."
    )

    parser.add_argument(
        "--optimize",
        action="store_true",
        help="Run module 6 dose optimization instead of a single simulation."
    )
    parser.add_argument(
        "--scan_output_dir",
        type=str,
        default="outputs/run_module6",
        help="Directory for dose-scan CSV and optimization figures."
    )

    return parser.parse_args()


def run_single_simulation(
    project_root: Path,
    config: dict,
    dose: float,
    t_end: float,
    dt: float,
    output_dir: str,
) -> None:
    results = run_simulation(
        config=config,
        dose=dose,
        t_end=t_end,
        dt=dt,
    )

    out_dir = ensure_dir(project_root / output_dir)
    figure_paths = plot_all(results, out_dir)
    save_results(results, out_dir)

    print("Single simulation finished.")
    print(f"Route: {config['route']}")
    print(f"Dose: {dose}")
    print(f"t_end: {t_end}")
    print(f"dt: {dt}")
    print(f"Output directory: {out_dir}")

    if not results.get("success", True):
        print(f"WARNING: solver reported failure: {results.get('message', 'Unknown error')}")

    print("\nGenerated figures:")
    for name, path in figure_paths.items():
        print(f"  {name}: {path}")


def run_optimization(
    project_root: Path,
    config: dict,
    t_end: float,
    dt: float,
    scan_output_dir: str,
) -> None:
    from optimization.min_effective_dose import run_dose_scan, find_min_effective_dose
    from optimization.report_module6 import save_dose_scan_table, plot_dose_response

    scan_results = run_dose_scan(
        base_config=config,
        t_end=t_end,
        dt=dt,
    )
    best = find_min_effective_dose(scan_results)

    scan_out_dir = ensure_dir(project_root / scan_output_dir)
    csv_path = save_dose_scan_table(scan_results, scan_out_dir)
    fig_paths = plot_dose_response(scan_results, scan_out_dir)

    print("Dose optimization finished.")
    print(f"Route: {config['route']}")
    print(f"t_end: {t_end}")
    print(f"dt: {dt}")
    print(f"Scan output directory: {scan_out_dir}")
    print(f"Dose scan saved to: {csv_path}")

    if best is None:
        print("\nNo feasible dose found under current thresholds.")
    else:
        print("\nMinimum feasible dose found:")
        for k, v in best.items():
            print(f"  {k}: {v}")

    print("\nGenerated optimization figures:")
    for name, path in fig_paths.items():
        print(f"  {name}: {path}")


def main() -> None:
    args = parse_args()
    project_root = Path(__file__).resolve().parent

    # Load base configuration
    config = load_base_config(project_root)

    # Override route if provided
    if args.route is not None:
        config["route"] = args.route

    # Resolve simulation settings EARLY so both branches can use them
    t_end = args.t_end if args.t_end is not None else config["simulation"]["default_t_end"]
    dt = args.dt if args.dt is not None else config["simulation"]["default_dt"]

    if args.optimize:
        run_optimization(
            project_root=project_root,
            config=config,
            t_end=t_end,
            dt=dt,
            scan_output_dir=args.scan_output_dir,
        )
    else:
        run_single_simulation(
            project_root=project_root,
            config=config,
            dose=args.dose,
            t_end=t_end,
            dt=dt,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()