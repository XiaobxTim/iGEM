from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

import numpy as np

from models.editing.module5 import compute_editing_metrics
from models.editing.offtarget_panel import (
    apply_offtarget_rows_to_config,
    summarize_offtarget_panel,
)
from models.editing.sequence_to_kinetics import (
    apply_design_to_config,
    load_design_table,
)
from models.full_model.simulator import run_simulation
from optimization.metrics import extract_summary_metrics
from optimization.min_effective_dose import find_min_effective_dose, run_dose_scan
from utils.config_loader import load_base_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DESIGN_TABLE = PROJECT_ROOT / "wetlab" / "templates" / "puf_design_candidates_template.csv"
MAX_PANEL_BYTES = 5 * 1024 * 1024
MAX_PANEL_ROWS = 10_000
PANEL_REQUIRED_COLUMNS = (
    "site_id",
    "gene",
    "initial_pool",
    "binding_score",
    "accessibility",
    "context_score",
    "validation_priority",
)
SCORE_COLUMNS = ("initial_pool", "binding_score", "accessibility", "context_score")
ALLOWED_MODES = {"single", "optimize"}
ALLOWED_ROUTES = {"footpad", "im", "iv"}
ALLOWED_DURATIONS = {24.0, 48.0, 72.0, 168.0}


def get_designs() -> list[dict[str, float | str]]:
    return load_design_table(DESIGN_TABLE)


def parse_candidate_panel(content: bytes) -> list[dict[str, float | str]]:
    """Validate and parse a PUF candidate-panel CSV without persisting uploads."""
    if len(content) > MAX_PANEL_BYTES:
        raise ValueError("Candidate panel exceeds the 5 MB upload limit.")
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("Candidate panel must be UTF-8 CSV.") from error

    reader = csv.DictReader(io.StringIO(text))
    fields = set(reader.fieldnames or ())
    missing = set(PANEL_REQUIRED_COLUMNS) - fields
    if missing:
        raise ValueError(f"Candidate panel is missing columns: {', '.join(sorted(missing))}.")

    rows: list[dict[str, float | str]] = []
    site_ids: set[str] = set()
    for line_number, raw in enumerate(reader, start=2):
        if len(rows) >= MAX_PANEL_ROWS:
            raise ValueError(f"Candidate panel exceeds the {MAX_PANEL_ROWS:,}-row limit.")
        site_id = str(raw.get("site_id", "")).strip()
        if not site_id:
            raise ValueError(f"site_id is required on CSV line {line_number}.")
        if site_id in site_ids:
            raise ValueError(f"Duplicate site_id '{site_id}' on CSV line {line_number}.")
        site_ids.add(site_id)
        row: dict[str, float | str] = {
            "site_id": site_id,
            "gene": str(raw.get("gene", "")).strip(),
            "validation_priority": str(raw.get("validation_priority", "")).strip(),
        }
        for column in SCORE_COLUMNS:
            try:
                value = float(str(raw.get(column, "")))
            except ValueError as error:
                raise ValueError(
                    f"{column} must be numeric on CSV line {line_number}."
                ) from error
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(
                    f"{column} must be between 0 and 1 on CSV line {line_number}."
                )
            row[column] = value
        rows.append(row)
    return rows


def _native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def validate_inputs(
    *,
    mode: str,
    design_id: str,
    route: str,
    dose: float,
    duration: float,
) -> None:
    if not 0.1 <= dose <= 10.0:
        raise ValueError("Dose must be between 0.1 and 10.0 normalized units.")
    if mode not in ALLOWED_MODES:
        raise ValueError("Mode must be 'single' or 'optimize'.")
    if route not in ALLOWED_ROUTES:
        raise ValueError("Route must be footpad, im, or iv.")
    if duration not in ALLOWED_DURATIONS:
        raise ValueError("Duration must be one of 24, 48, 72, or 168 hours.")
    if design_id not in {str(row["design_id"]) for row in get_designs()}:
        raise ValueError(f"Unknown design preset '{design_id}'.")


def _prepare_config(
    design_id: str,
    route: str,
    panel_rows: list[dict[str, float | str]] | None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    design = next(row for row in get_designs() if row["design_id"] == design_id)
    config = load_base_config(PROJECT_ROOT)
    config["route"] = route
    config, modifiers = apply_design_to_config(config, design)
    if panel_rows is None:
        panel_summary: dict[str, Any] = {
            "provided": False,
            "n_sites": 0.0,
            "interpretation": "Base-config aggregate off-target prior retained.",
        }
    elif not panel_rows:
        panel_summary = {
            "provided": True,
            "n_sites": 0.0,
            "warning": (
                "The uploaded panel contained no candidate rows; the conservative "
                "base-config aggregate off-target prior was retained."
            ),
        }
    else:
        config, panel_summary = apply_offtarget_rows_to_config(config, panel_rows)
        panel_summary["provided"] = True
    return config, dict(modifiers), panel_summary


def run_model(
    *,
    mode: str,
    design_id: str,
    route: str,
    dose: float,
    duration: float,
    panel_rows: list[dict[str, float | str]] | None = None,
) -> dict[str, Any]:
    validate_inputs(
        mode=mode,
        design_id=design_id,
        route=route,
        dose=dose,
        duration=duration,
    )
    config, modifiers, panel_summary = _prepare_config(design_id, route, panel_rows)
    common = {
        "mode": mode,
        "inputs": {
            "design_id": design_id,
            "route": route,
            "dose": dose,
            "duration_hours": duration,
        },
        "design_modifiers": modifiers,
        "panel_summary": panel_summary,
        "disclaimer": (
            "Literature-informed mechanistic simulation; not clinically calibrated and not "
            "a prediction of patient outcome. Experimental validation is required."
        ),
    }

    if mode == "optimize":
        scan = run_dose_scan(config, t_end=duration, dt=0.2)
        return _native(
            {
                **common,
                "dose_scan": scan,
                "minimum_feasible_dose": find_min_effective_dose(scan),
            }
        )

    result = run_simulation(config, dose=dose, t_end=duration, dt=0.2)
    if not result["success"]:
        raise RuntimeError(f"Numerical solver failed: {result['message']}")
    metrics = compute_editing_metrics(result["y"], result["idx"], config)
    summary = extract_summary_metrics(result)
    series = {
        "time": result["t"],
        "pbrain": result["y"][result["idx"]["P_brain"]],
        "apoe3_like": metrics["apoe3_like_fraction"],
        "apoe2_like": metrics["apoe2_like_fraction"],
        "off_target": metrics["off_target_burden"],
        "specificity": metrics["specificity_index"],
    }
    return _native({**common, "metrics": summary, "series": series})


def baseline_panel_summary() -> dict[str, float]:
    """Expose the empty-panel aggregate for documentation and tests."""
    return summarize_offtarget_panel([])
