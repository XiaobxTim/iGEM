from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List
import csv
import math

import numpy as np


@dataclass(frozen=True)
class WetlabObservation:
    """
    One measured wet-lab data point.

    The model itself works in arbitrary normalized units because the current
    project is still pre-calibration. This bridge lets us compare simulated
    trajectories to assay-level observations once the wet lab starts producing
    qPCR, amplicon sequencing, Western blot, or viability data.
    """

    assay: str
    time_h: float
    value: float
    sd: float = 1.0
    weight: float = 1.0


def _interp_state(results: Dict, state_name: str, time_h: float) -> float:
    """Interpolate one model state at the assay sampling time."""
    t = np.asarray(results["t"], dtype=float)
    y = np.asarray(results["y"], dtype=float)
    idx = results["idx"]
    return float(np.interp(float(time_h), t, y[idx[state_name], :]))


def simulated_assay_value(results: Dict, assay: str, time_h: float) -> float:
    """
    Convert model states into wet-lab observable quantities.

    Supported assays intentionally match realistic experiments a student team
    can run:
    - vector_qpcr_blood: vector genomes in blood by qPCR/ddPCR
    - vector_qpcr_liver: liver vector burden proxy
    - vector_qpcr_brain: brain/nuclear vector delivery proxy
    - puf_mrna_qpcr: expressed PUF-APOBEC mRNA
    - editor_western: active editor protein proxy
    - on_editing_amplicon_pct: on-target editing percent by amplicon sequencing
    - off_editing_amplicon_pct: off-target editing percent by amplicon sequencing
    - viability_pct: simple burden proxy, to be replaced by measured viability
    """
    eps = 1e-12
    assay = assay.strip()

    if assay == "vector_qpcr_blood":
        return _interp_state(results, "A_blood", time_h)
    if assay == "vector_qpcr_liver":
        return _interp_state(results, "A_liver", time_h)
    if assay == "vector_qpcr_brain":
        return _interp_state(results, "A_brain_nuc", time_h)
    if assay == "puf_mrna_qpcr":
        return _interp_state(results, "mRNA_brain", time_h)
    if assay == "editor_western":
        return _interp_state(results, "P_brain", time_h)

    if assay == "on_editing_amplicon_pct":
        if "S_APOE4" in results["idx"]:
            apoe4 = _interp_state(results, "S_APOE4", time_h)
            apoe3 = _interp_state(results, "S_APOE3_like", time_h)
            apoe2 = _interp_state(results, "S_APOE2_like", time_h)
            mixed = _interp_state(results, "S_APOE158_only", time_h)
            total = apoe4 + apoe3 + apoe2 + mixed + eps
            return 100.0 * (apoe3 + apoe2) / total
        s_on = _interp_state(results, "S_on", time_h)
        e_on = _interp_state(results, "E_on", time_h)
        return 100.0 * e_on / (s_on + e_on + eps)

    if assay == "apoe3_like_amplicon_pct":
        apoe4 = _interp_state(results, "S_APOE4", time_h)
        apoe3 = _interp_state(results, "S_APOE3_like", time_h)
        apoe2 = _interp_state(results, "S_APOE2_like", time_h)
        mixed = _interp_state(results, "S_APOE158_only", time_h)
        return 100.0 * apoe3 / (apoe4 + apoe3 + apoe2 + mixed + eps)

    if assay == "apoe2_like_amplicon_pct":
        apoe4 = _interp_state(results, "S_APOE4", time_h)
        apoe3 = _interp_state(results, "S_APOE3_like", time_h)
        apoe2 = _interp_state(results, "S_APOE2_like", time_h)
        mixed = _interp_state(results, "S_APOE158_only", time_h)
        return 100.0 * apoe2 / (apoe4 + apoe3 + apoe2 + mixed + eps)

    if assay == "off_editing_amplicon_pct":
        if "B_local_bystander" in results["idx"]:
            local = _interp_state(results, "B_local_bystander", time_h)
            puf_off = _interp_state(results, "E_puf_off", time_h)
            deam = _interp_state(results, "E_deaminase_bg", time_h)
            apoe4 = _interp_state(results, "S_APOE4", time_h)
            apoe3 = _interp_state(results, "S_APOE3_like", time_h)
            apoe2 = _interp_state(results, "S_APOE2_like", time_h)
            mixed = _interp_state(results, "S_APOE158_only", time_h)
            total = apoe4 + apoe3 + apoe2 + mixed + local + puf_off + deam + eps
            return 100.0 * (local + puf_off + deam) / total
        s_off = _interp_state(results, "S_off", time_h)
        e_off = _interp_state(results, "E_off", time_h)
        return 100.0 * e_off / (s_off + e_off + eps)

    if assay == "local_bystander_amplicon_pct":
        local = _interp_state(results, "B_local_bystander", time_h)
        apoe4 = _interp_state(results, "S_APOE4", time_h)
        apoe3 = _interp_state(results, "S_APOE3_like", time_h)
        apoe2 = _interp_state(results, "S_APOE2_like", time_h)
        mixed = _interp_state(results, "S_APOE158_only", time_h)
        return 100.0 * local / (apoe4 + apoe3 + apoe2 + mixed + local + eps)

    if assay == "puf_mismatch_offtarget_pct":
        s_off = _interp_state(results, "S_puf_off", time_h)
        e_off = _interp_state(results, "E_puf_off", time_h)
        return 100.0 * e_off / (s_off + e_off + eps)

    if assay == "deaminase_background_pct":
        s_bg = _interp_state(results, "S_deaminase_bg", time_h)
        e_bg = _interp_state(results, "E_deaminase_bg", time_h)
        return 100.0 * e_bg / (s_bg + e_bg + eps)

    if assay == "viability_pct":
        # A first-pass, interpretable burden proxy. Liver exposure and high
        # editor expression lower the score. Wet-lab viability data can replace
        # this proxy by fitting the coefficients later.
        liver = _interp_state(results, "A_liver", time_h)
        editor = _interp_state(results, "P_brain", time_h)
        burden = 0.15 * liver + 0.05 * editor
        return 100.0 * math.exp(-burden)

    raise KeyError(f"Unknown wet-lab assay '{assay}'.")


def load_wetlab_observations(path: str | Path) -> List[WetlabObservation]:
    """Load observations from the CSV template in wetlab/templates/."""
    observations: List[WetlabObservation] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            observations.append(
                WetlabObservation(
                    assay=row["assay"],
                    time_h=float(row["time_h"]),
                    value=float(row["value"]),
                    sd=float(row.get("sd") or 1.0),
                    weight=float(row.get("weight") or 1.0),
                )
            )
    return observations


def residual_table(results: Dict, observations: Iterable[WetlabObservation]) -> List[Dict[str, float | str]]:
    """
    Build a row-by-row residual table for calibration and wiki reporting.
    """
    rows: List[Dict[str, float | str]] = []
    for obs in observations:
        pred = simulated_assay_value(results, obs.assay, obs.time_h)
        sigma = max(float(obs.sd), 1e-12)
        residual = (pred - obs.value) / sigma
        rows.append(
            {
                "assay": obs.assay,
                "time_h": obs.time_h,
                "observed": obs.value,
                "predicted": pred,
                "sd": obs.sd,
                "weight": obs.weight,
                "weighted_residual": obs.weight * residual,
            }
        )
    return rows


def weighted_rmse(results: Dict, observations: Iterable[WetlabObservation]) -> float:
    """
    A calibration objective that can drive future parameter fitting.

    Lower values mean the simulated delivery-expression-editing chain is closer
    to the wet-lab measurements after accounting for assay uncertainty.
    """
    rows = residual_table(results, observations)
    if not rows:
        return float("nan")
    squared = [float(row["weighted_residual"]) ** 2 for row in rows]
    return float(math.sqrt(sum(squared) / len(squared)))
