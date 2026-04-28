from pathlib import Path

from utils.config_loader import load_base_config
from optimization.min_effective_dose import run_dose_scan, find_min_effective_dose
from optimization.report_module6 import save_dose_scan_table, plot_dose_response


def main():
    project_root = Path(__file__).resolve().parent
    config = load_base_config(project_root)

    scan_results = run_dose_scan(base_config=config)
    best = find_min_effective_dose(scan_results)

    output_dir = project_root / "outputs" / "run_module6"
    csv_path = save_dose_scan_table(scan_results, output_dir)
    fig_paths = plot_dose_response(scan_results, output_dir)

    print(f"Dose scan saved to: {csv_path}")
    if best is None:
        print("No feasible dose found under current thresholds.")
    else:
        print("Minimum feasible dose found:")
        for k, v in best.items():
            print(f"  {k}: {v}")

    for name, path in fig_paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()