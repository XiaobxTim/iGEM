from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import numpy as np

from spatial_model.intracellular_reference import simulate_intracellular
from spatial_model.parameter_contract import build_spatial_parameters, write_parameter_file
from utils.config_loader import load_base_config


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CPP_CORE = PROJECT_ROOT / "spatial_model" / "cpp_core"


def test_cpp_rk4_matches_python_reference(tmp_path):
    build = subprocess.run(
        ["make", "-C", str(CPP_CORE), "clean", "all"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert build.returncode == 0, build.stdout + build.stderr

    params = build_spatial_parameters(load_base_config("."))
    parameter_file = write_parameter_file(params, tmp_path / "parameters.cfg")
    run = subprocess.run(
        [
            str(CPP_CORE / "ode_probe"),
            "--parameters",
            str(parameter_file),
            "--minutes",
            "120",
            "--dt",
            "1",
            "--uptake",
            "0.0005",
            "--apoe-scale",
            "1",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    rows = list(csv.DictReader(run.stdout.splitlines()))
    cpp_final = {key: float(value) for key, value in rows[-1].items() if key != "time_min"}

    python = simulate_intracellular(
        params,
        minutes=120.0,
        dt_min=1.0,
        uptake_amount_per_min=0.0005,
        apoe_scale=1.0,
    )
    assert set(cpp_final) == set(python.state_names)
    reference_final = python.values[-1]
    cpp_values = np.array([cpp_final[name] for name in python.state_names])
    scale = np.maximum(np.abs(reference_final), 1e-8)
    assert np.max(np.abs(cpp_values - reference_final) / scale) < 1e-3


def test_cpp_probe_rejects_negative_input(tmp_path):
    subprocess.run(
        ["make", "-C", str(CPP_CORE), "all"],
        text=True,
        capture_output=True,
        check=True,
    )
    parameter_file = write_parameter_file(
        build_spatial_parameters(load_base_config(".")),
        tmp_path / "parameters.cfg",
    )

    run = subprocess.run(
        [
            str(CPP_CORE / "ode_probe"),
            "--parameters",
            str(parameter_file),
            "--uptake",
            "-1",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert run.returncode != 0
    assert "uptake must be non-negative" in run.stderr
