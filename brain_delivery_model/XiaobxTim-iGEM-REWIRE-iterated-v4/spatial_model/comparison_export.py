from __future__ import annotations

import argparse
from pathlib import Path

from spatial_model.paraview_export import write_comparison_series


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a side-by-side ParaView dose series")
    parser.add_argument(
        "--scenario",
        action="append",
        required=True,
        help="Scenario in LABEL=/absolute/path/to/simulation.pvd form",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    arguments = parser.parse_args()
    scenarios: list[tuple[str, Path]] = []
    for value in arguments.scenario:
        if "=" not in value:
            parser.error("--scenario requires LABEL=PVD")
        label, filename = value.split("=", 1)
        scenarios.append((label, Path(filename)))
    result = write_comparison_series(scenarios, arguments.output_dir)
    print(f"Comparison series written to {result}")


if __name__ == "__main__":
    main()
