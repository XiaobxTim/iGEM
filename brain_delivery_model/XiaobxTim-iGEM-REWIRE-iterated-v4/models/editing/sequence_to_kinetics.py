from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, List
import csv
import math


def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def load_design_table(path: str | Path) -> List[Dict[str, float | str]]:
    """
    Load PUF-deaminase design candidates.

    The table is intentionally lightweight: it can be filled from hand-curated
    APOE target designs now and replaced by REWIRE/ProAPOBEC benchmark data
    later.
    """
    rows: List[Dict[str, float | str]] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            rows.append(
                {
                    "design_id": raw["design_id"],
                    "editor_type": raw.get("editor_type", "A3A"),
                    "puf_repeats": float(raw.get("puf_repeats") or 10),
                    "target_site": raw.get("target_site", "APOE112"),
                    "puf_score_112": float(raw.get("puf_score_112") or 0.0),
                    "puf_score_158": float(raw.get("puf_score_158") or 0.0),
                    "mismatch_count": float(raw.get("mismatch_count") or 0.0),
                    "accessibility_112": float(raw.get("accessibility_112") or 0.5),
                    "accessibility_158": float(raw.get("accessibility_158") or 0.5),
                    "distance_112": float(raw.get("distance_112") or 2.0),
                    "distance_158": float(raw.get("distance_158") or 2.0),
                    "uc_context_112": float(raw.get("uc_context_112") or 1.0),
                    "uc_context_158": float(raw.get("uc_context_158") or 0.5),
                    "local_bystander_risk": float(raw.get("local_bystander_risk") or 0.05),
                    "deaminase_background_prior": float(raw.get("deaminase_background_prior") or 1.0),
                    "editor_activity_scale": float(raw.get("editor_activity_scale") or 1.0),
                }
            )
    return rows


def design_to_kinetic_modifiers(row: Dict[str, float | str]) -> Dict[str, float | str]:
    """
    Map sequence/design features to kinetic modifiers.

    This implements the optimization-plan idea:
    sequence/accessibility -> k_on, context/distance -> k_cat.
    """
    puf_repeats = float(row["puf_repeats"])
    mismatch_count = float(row["mismatch_count"])
    repeat_bonus = 1.0 + 0.10 * max(puf_repeats - 8.0, 0.0)
    mismatch_penalty = math.exp(-0.35 * mismatch_count)
    # Han et al. NAR 2022 reported a strong specificity gain when moving from
    # 8-repeat to 10-repeat PUF recognition. We encode that as a separate
    # off-target-risk modifier instead of only increasing on-target binding.
    puf_offtarget_scale = math.exp(-1.4 * max(puf_repeats - 8.0, 0.0)) * math.exp(
        0.45 * mismatch_count
    )

    kon112_scale = math.exp(
        0.9 * float(row["puf_score_112"])
        + 0.7 * float(row["accessibility_112"])
    ) * repeat_bonus * mismatch_penalty
    kon158_scale = math.exp(
        0.9 * float(row["puf_score_158"])
        + 0.7 * float(row["accessibility_158"])
    ) * repeat_bonus * mismatch_penalty

    distance112 = abs(float(row["distance_112"]) - 2.0)
    distance158 = abs(float(row["distance_158"]) - 2.0)
    # CU-REWIRE has a narrow C-to-U window, strongest around the second
    # nucleotide after the PUF binding site, and a strong UC-context preference.
    # The template stores context as a 0-1 relative motif score.
    kcat112_scale = (0.04 + 0.96 * float(row["uc_context_112"])) * math.exp(
        -0.9 * distance112
    )
    kcat158_scale = (0.04 + 0.96 * float(row["uc_context_158"])) * math.exp(
        -0.9 * distance158
    )
    editor_activity_scale = max(float(row.get("editor_activity_scale", 1.0)), 0.0)

    return {
        "design_id": row["design_id"],
        "editor_type": row["editor_type"],
        "kon112_scale": kon112_scale,
        "kon158_scale": kon158_scale,
        "kcat112_scale": kcat112_scale * editor_activity_scale,
        "kcat158_scale": kcat158_scale * editor_activity_scale,
        "local_bystander_scale": max(float(row["local_bystander_risk"]), 0.0),
        "deaminase_background_scale": max(float(row["deaminase_background_prior"]), 0.0)
        * editor_activity_scale,
        "editor_activity_scale": editor_activity_scale,
        "puf_offtarget_scale": max(puf_offtarget_scale, 0.0),
    }


def apply_design_to_config(config: Dict, row: Dict[str, float | str]) -> tuple[Dict, Dict[str, float | str]]:
    updated = deepcopy(config)
    modifiers = design_to_kinetic_modifiers(row)
    editing = updated["editing"]

    editing["editor_type"] = str(modifiers["editor_type"])
    editing["k_on_112"] = float(editing["k_on_112"]) * float(modifiers["kon112_scale"])
    editing["k_on_158"] = float(editing["k_on_158"]) * float(modifiers["kon158_scale"])
    editing["k_cat_112"] = float(editing["k_cat_112"]) * float(modifiers["kcat112_scale"])
    editing["k_cat_158"] = float(editing["k_cat_158"]) * float(modifiers["kcat158_scale"])
    editing["local_bystander_per_112"] = float(modifiers["local_bystander_scale"])
    editing["local_bystander_per_158"] = 1.5 * float(modifiers["local_bystander_scale"])
    editing["k_on_puf_off"] = float(editing["k_on_puf_off"]) * float(
        modifiers["puf_offtarget_scale"]
    )
    editing["k_cat_puf_off"] = float(editing["k_cat_puf_off"]) * float(
        modifiers["puf_offtarget_scale"]
    )
    editing["k_deaminase_bg"] = float(editing["k_deaminase_bg"]) * float(
        modifiers["deaminase_background_scale"]
    )
    return updated, modifiers


def pareto_front(rows: Iterable[Dict[str, float]]) -> List[Dict[str, float]]:
    """
    Return nondominated design rows.

    We maximize utility, APOE3-like fraction and specificity; minimize
    off-target burden and LDLR-risk proxy.
    """
    rows = list(rows)
    front = []
    for row in rows:
        dominated = False
        for other in rows:
            if row is other:
                continue
            no_worse = (
                other["utility"] >= row["utility"]
                and other["apoe3_like_fraction_final"] >= row["apoe3_like_fraction_final"]
                and other["specificity_index_final"] >= row["specificity_index_final"]
                and other["off_target_burden_final"] <= row["off_target_burden_final"]
                and other["ldlr_binding_risk_proxy_final"] <= row["ldlr_binding_risk_proxy_final"]
            )
            strictly_better = (
                other["utility"] > row["utility"]
                or other["apoe3_like_fraction_final"] > row["apoe3_like_fraction_final"]
                or other["specificity_index_final"] > row["specificity_index_final"]
                or other["off_target_burden_final"] < row["off_target_burden_final"]
                or other["ldlr_binding_risk_proxy_final"] < row["ldlr_binding_risk_proxy_final"]
            )
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            front.append(row)
    return front
