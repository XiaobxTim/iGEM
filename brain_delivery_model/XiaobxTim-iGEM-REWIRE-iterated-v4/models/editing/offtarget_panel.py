from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, List
import csv


def load_offtarget_panel(path: str | Path) -> List[Dict[str, float | str]]:
    """
    Load a REWIRE-style off-target candidate panel.

    Each row represents one candidate RNA site. The model does not yet track
    every site as a separate ODE state, so this module compresses the panel into
    an effective aggregate off-target pool. That is a pragmatic bridge between
    transcriptome-wide dry-lab screening and the current compact ODE system.
    """
    rows: List[Dict[str, float | str]] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            rows.append(
                {
                    "site_id": raw["site_id"],
                    "gene": raw.get("gene", ""),
                    "initial_pool": float(raw.get("initial_pool") or 1.0),
                    "binding_score": float(raw.get("binding_score") or 0.0),
                    "accessibility": float(raw.get("accessibility") or 1.0),
                    "context_score": float(raw.get("context_score") or 1.0),
                    "validation_priority": raw.get("validation_priority", "medium"),
                }
            )
    return rows


def summarize_offtarget_panel(rows: Iterable[Dict[str, float | str]]) -> Dict[str, float]:
    """
    Convert candidate-site rows to effective aggregate model quantities.

    Interpretation:
    - `effective_pool` becomes `editing.S_off_init`.
    - `relative_binding` scales `editing.k_on_off`.
    - `relative_catalysis` scales `editing.k_cat_off`.

    The weights intentionally use accessibility because hidden RNA sites should
    contribute less to binding and editing burden.
    """
    rows = list(rows)
    if not rows:
        return {
            "n_sites": 0,
            "effective_pool": 0.0,
            "relative_binding": 0.0,
            "relative_catalysis": 0.0,
        }

    weighted_pool = 0.0
    weighted_binding = 0.0
    weighted_catalysis = 0.0

    for row in rows:
        pool = float(row["initial_pool"])
        accessibility = max(float(row["accessibility"]), 0.0)
        binding_score = max(float(row["binding_score"]), 0.0)
        context_score = max(float(row["context_score"]), 0.0)
        weight = pool * accessibility

        weighted_pool += weight
        weighted_binding += weight * binding_score
        weighted_catalysis += weight * binding_score * context_score

    if weighted_pool <= 0:
        relative_binding = 0.0
        relative_catalysis = 0.0
    else:
        relative_binding = weighted_binding / weighted_pool
        relative_catalysis = weighted_catalysis / weighted_pool

    return {
        "n_sites": float(len(rows)),
        "effective_pool": weighted_pool,
        "relative_binding": relative_binding,
        "relative_catalysis": relative_catalysis,
    }


def apply_offtarget_panel_to_config(
    config: Dict,
    panel_path: str | Path,
    min_pool: float = 1e-9,
) -> tuple[Dict, Dict[str, float]]:
    """
    Return a config copy with aggregate off-target parameters updated.

    This keeps the current ODE model compact while letting the off-target risk
    come from a transparent candidate table instead of a hand-picked `S_off`.
    """
    updated = deepcopy(config)
    rows = load_offtarget_panel(panel_path)
    summary = summarize_offtarget_panel(rows)
    editing = updated["editing"]

    if editing.get("model") == "apoe_multisite":
        pool_key = "S_puf_off_init"
        kon_key = "k_on_puf_off"
        kcat_key = "k_cat_puf_off"
    else:
        pool_key = "S_off_init"
        kon_key = "k_on_off"
        kcat_key = "k_cat_off"

    base_pool = max(float(editing.get(pool_key, 1.0)), min_pool)
    base_k_on = float(editing.get(kon_key, 0.0))
    base_k_cat = float(editing.get(kcat_key, 0.0))

    effective_pool = max(summary["effective_pool"], min_pool)
    editing[pool_key] = effective_pool

    # Keep the original parameter values as priors, then scale by panel-derived
    # relative risk. The clamp avoids accidentally erasing off-target risk when
    # an early placeholder panel has sparse scores.
    editing[kon_key] = base_k_on * max(summary["relative_binding"], 0.05)
    editing[kcat_key] = base_k_cat * max(summary["relative_catalysis"], 0.05)

    summary["pool_key"] = pool_key
    summary["base_offtarget_pool"] = base_pool
    summary["applied_offtarget_pool"] = editing[pool_key]
    summary["applied_offtarget_k_on"] = editing[kon_key]
    summary["applied_offtarget_k_cat"] = editing[kcat_key]

    return updated, summary
