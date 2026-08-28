from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from spatial_model.parameter_contract import build_spatial_parameters, write_parameter_file
from spatial_model.paraview_export import publish_atomically
from spatial_model.pbpk_boundary import build_boundary_curves, write_boundary_csv
from spatial_model.physicell_stage import stage_physicell_project
from spatial_model.project_builder import render_physicell_config
from spatial_model.raw_output import discover_snapshots, load_snapshot
from spatial_model.tissue_geometry import generate_tissue, write_cells_csv
from utils.config_loader import load_base_config


DEFAULT_DOSES = "0.3,1.0,3.0"
DEFAULT_PHYSICELL = Path("/private/tmp/PhysiCell-1.14.2/PhysiCell")
DEFAULT_PVPYTHON = Path("/Applications/ParaView-6.0.1.app/Contents/bin/pvpython")


def parse_doses(value: str) -> tuple[float, float, float]:
    try:
        doses = tuple(float(part.strip()) for part in value.split(",") if part.strip())
    except ValueError as error:
        raise ValueError("doses must be three comma-separated numbers") from error
    if (
        len(doses) != 3
        or len(set(doses)) != 3
        or any(not math.isfinite(dose) or dose <= 0.0 for dose in doses)
        or 1.0 not in doses
    ):
        raise ValueError("doses must be three unique positive values including 1.0")
    return doses  # type: ignore[return-value]


def scenario_name(dose: float) -> str:
    return f"dose_{dose:.1f}".replace(".", "p")


def _run(command: list[str], *, cwd: Path, log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        tail = "\n".join(log.read_text(encoding="utf-8").splitlines()[-30:])
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{tail}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _collect_metrics(staging: Path, scenarios: list[tuple[float, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for dose, name in scenarios:
        metrics_path = staging / name / "paraview" / "metrics.csv"
        with metrics_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                rows.append({"scenario": name, "dose": dose, **row})
    output = staging / "dose_comparison_metrics.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _validate_mass_balance(raw_dir: Path) -> float:
    with (raw_dir / "mass_balance.csv").open(newline="", encoding="utf-8") as handle:
        errors = [float(row["relative_error"]) for row in csv.DictReader(handle)]
    maximum = max(errors, default=math.inf)
    if not math.isfinite(maximum) or maximum > 0.01:
        raise RuntimeError(f"mass-balance relative error exceeds 1%: {maximum:.6g}")
    return maximum


def _scientific_checks(
    staging: Path, scenarios: list[tuple[float, str]], metric_rows: list[dict[str, object]]
) -> dict[str, object]:
    final_metrics: dict[str, dict[str, float]] = {}
    spatial: dict[str, dict[str, float | bool]] = {}
    for _, name in scenarios:
        rows = [row for row in metric_rows if row["scenario"] == name]
        final_metrics[name] = {key: float(value) for key, value in rows[-1].items() if key not in {"scenario"}}
        snapshot = load_snapshot(discover_snapshots(staging / name / "raw")[-1])
        types = snapshot.cell_data["cell_type"]
        editing = snapshot.cell_data["editing_fraction"]
        distance = snapshot.cell_data["distance_to_vessel_um"]
        brain = (types == 2) | (types == 3)
        low_quantile, high_quantile = np.quantile(distance[brain], [0.25, 0.75])
        near = brain & (distance <= low_quantile)
        far = brain & (distance >= high_quantile)
        neuron = types == 2
        astrocyte = types == 3
        near_mean = float(np.mean(editing[near]))
        far_mean = float(np.mean(editing[far]))
        neuron_mean = float(np.mean(editing[neuron]))
        astrocyte_mean = float(np.mean(editing[astrocyte]))
        spatial[name] = {
            "near_vessel_editing_mean": near_mean,
            "far_vessel_editing_mean": far_mean,
            "astrocyte_editing_mean": astrocyte_mean,
            "neuron_editing_mean": neuron_mean,
            "near_exceeds_far": near_mean >= far_mean,
            "astrocyte_exceeds_neuron": astrocyte_mean >= neuron_mean,
        }
    ordered = [name for _, name in sorted(scenarios)]
    mass = [final_metrics[name]["extracellular_AAV_mass"] for name in ordered]
    editing = [final_metrics[name]["brain_cell_editing_mean"] for name in ordered]
    dose_monotonic = all(right >= left for left, right in zip(mass, mass[1:])) and all(
        right >= left for left, right in zip(editing, editing[1:])
    )
    spatial_valid = all(
        values["near_exceeds_far"] and values["astrocyte_exceeds_neuron"]
        for values in spatial.values()
    )
    result = {
        "dose_response_monotonic": dose_monotonic,
        "spatial_priors_satisfied": spatial_valid,
        "final_spatial_metrics": spatial,
    }
    (staging / "scientific_checks.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    if not dose_monotonic or not spatial_valid:
        raise RuntimeError(f"scientific validation failed: {result}")
    return result


def run_pipeline(
    *,
    project_root: Path,
    physicell_root: Path,
    pvpython: Path,
    output_dir: Path,
    doses: tuple[float, float, float],
    duration_hours: float,
    boundary_dt_min: float,
    output_interval_min: float,
    jobs: int,
) -> Path:
    """Build, run, convert, validate, and atomically publish three scenarios."""

    if output_dir.exists():
        raise FileExistsError(output_dir)
    if not physicell_root.is_dir():
        raise FileNotFoundError(physicell_root)
    if not pvpython.is_file():
        raise FileNotFoundError(pvpython)
    if min(duration_hours, boundary_dt_min, output_interval_min) <= 0.0:
        raise ValueError("duration and output intervals must be positive")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}-", dir=output_dir.parent))
    try:
        config = load_base_config(project_root)
        parameters = build_spatial_parameters(config)
        curves = build_boundary_curves(
            config,
            doses,
            t_end_h=duration_hours,
            dt_h=boundary_dt_min / 60.0,
            reference_dose=1.0,
        )
        cells = generate_tissue(seed=42)

        stage_physicell_project(physicell_root, project_root / "spatial_model")
        _run(
            ["make", f"-j{max(1, jobs)}", "all"],
            cwd=physicell_root,
            log=staging / "logs" / "build.log",
        )
        executable = physicell_root / "brain_delivery"
        template = physicell_root / "sample_projects/template/config/PhysiCell_settings.xml"
        scenario_specs: list[tuple[float, str]] = []
        mass_errors: dict[str, float] = {}
        for dose in doses:
            name = scenario_name(dose)
            scenario_specs.append((dose, name))
            scenario = staging / name
            inputs = scenario / "inputs"
            raw = scenario / "raw"
            inputs.mkdir(parents=True)
            raw.mkdir()
            boundary = write_boundary_csv(curves[dose], inputs / "brain_blood_input.csv")
            parameter_file = write_parameter_file(parameters, inputs / "spatial_parameters.cfg")
            cells_file = write_cells_csv(cells, inputs / "cells.csv")
            render_physicell_config(
                template,
                scenario / "PhysiCell_settings.xml",
                boundary_csv=Path("inputs/brain_blood_input.csv"),
                parameter_file=Path("inputs/spatial_parameters.cfg"),
                cells_csv=Path("inputs/cells.csv"),
                output_dir=Path("raw"),
                parameters=parameters,
                max_time_min=duration_hours * 60.0,
                output_interval_min=output_interval_min,
                seed=42,
            )
            _run(
                [str(executable), "PhysiCell_settings.xml"],
                cwd=scenario,
                log=staging / "logs" / f"{name}_physicell.log",
            )
            mass_errors[name] = _validate_mass_balance(raw)
            _run(
                [
                    str(pvpython),
                    "-m",
                    "spatial_model.paraview_export",
                    "--raw-dir",
                    str(raw),
                    "--output-dir",
                    str(scenario / "paraview"),
                ],
                cwd=project_root,
                log=staging / "logs" / f"{name}_paraview_export.log",
            )
            _run(
                [
                    str(pvpython),
                    "-m",
                    "spatial_model.validate_paraview",
                    str(scenario / "paraview/simulation.pvd"),
                ],
                cwd=project_root,
                log=staging / "logs" / f"{name}_paraview_validate.log",
            )
            shutil.copy2(raw / "mass_balance.csv", scenario / "mass_balance.csv")

        comparison_command = [
            str(pvpython),
            "-m",
            "spatial_model.comparison_export",
        ]
        for _, name in scenario_specs:
            comparison_command.extend(
                ["--scenario", f"{name}={staging / name / 'paraview/simulation.pvd'}"]
            )
        comparison_command.extend(["--output-dir", str(staging / "comparison")])
        _run(
            comparison_command,
            cwd=project_root,
            log=staging / "logs" / "comparison_export.log",
        )
        _run(
            [
                str(pvpython),
                "-m",
                "spatial_model.validate_comparison",
                str(staging / "comparison/comparison.pvd"),
            ],
            cwd=project_root,
            log=staging / "logs" / "comparison_validate.log",
        )
        metric_rows = _collect_metrics(staging, scenario_specs)
        scientific_checks = _scientific_checks(staging, scenario_specs, metric_rows)
        metadata = {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "model": "PhysiCell-BioFVM spatial brain delivery",
            "physicell_version": (physicell_root / "VERSION.txt").read_text(encoding="utf-8").strip(),
            "paraview_pvpython": str(pvpython),
            "doses": list(doses),
            "duration_hours": duration_hours,
            "boundary_dt_min": boundary_dt_min,
            "output_interval_min": output_interval_min,
            "seed": 42,
            "cell_counts": {"endothelial": 280, "neuron": 800, "astrocyte": 400},
            "mass_balance_max_relative_error": mass_errors,
            "scientific_checks": scientific_checks,
            "input_sha256": {
                name: {
                    "boundary": _sha256(staging / name / "inputs/brain_blood_input.csv"),
                    "parameters": _sha256(staging / name / "inputs/spatial_parameters.cfg"),
                    "cells": _sha256(staging / name / "inputs/cells.csv"),
                }
                for _, name in scenario_specs
            },
        }
        (staging / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        (staging / "OPEN_IN_PARAVIEW.txt").write_text(
            "Open comparison/comparison.pvd for the side-by-side dose comparison.\n"
            "Open dose_0p3/paraview/simulation.pvd, dose_1p0/paraview/simulation.pvd, "
            "or dose_3p0/paraview/simulation.pvd for an individual simulation.\n",
            encoding="utf-8",
        )
        return publish_atomically(staging, output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Run the PhysiCell brain-delivery pipeline")
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--physicell-root", type=Path, default=DEFAULT_PHYSICELL)
    parser.add_argument("--pvpython", type=Path, default=DEFAULT_PVPYTHON)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "outputs/physicell_brain_delivery",
    )
    parser.add_argument("--doses", default=DEFAULT_DOSES)
    parser.add_argument("--duration-hours", type=float, default=72.0)
    parser.add_argument("--boundary-dt-min", type=float, default=5.0)
    parser.add_argument("--output-interval-min", type=float, default=120.0)
    parser.add_argument("--jobs", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    arguments = _parse_args()
    result = run_pipeline(
        project_root=arguments.project_root.resolve(),
        physicell_root=arguments.physicell_root.resolve(),
        pvpython=arguments.pvpython.resolve(),
        output_dir=arguments.output_dir.resolve(),
        doses=parse_doses(arguments.doses),
        duration_hours=arguments.duration_hours,
        boundary_dt_min=arguments.boundary_dt_min,
        output_interval_min=arguments.output_interval_min,
        jobs=arguments.jobs,
    )
    print(f"Completed: {result}")
    print(f"Open in ParaView: {result / 'comparison/comparison.pvd'}")


if __name__ == "__main__":
    main()
