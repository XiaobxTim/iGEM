from __future__ import annotations

import numpy as np

from models.calibration.wetlab_bridge import simulated_assay_value


def _fake_results():
    state_order = [
        "A_blood",
        "A_liver",
        "A_brain_nuc",
        "mRNA_brain",
        "P_brain",
        "S_on",
        "E_on",
        "S_off",
        "E_off",
    ]
    idx = {name: i for i, name in enumerate(state_order)}
    y = np.array(
        [
            [1.0, 0.5],
            [0.0, 2.0],
            [0.0, 0.1],
            [0.0, 0.2],
            [0.0, 0.3],
            [1.0, 0.8],
            [0.0, 0.2],
            [5.0, 4.9],
            [0.0, 0.1],
        ]
    )
    return {"t": np.array([0.0, 10.0]), "y": y, "idx": idx}


def test_amplicon_percent_observables_are_state_derived():
    results = _fake_results()
    assert abs(simulated_assay_value(results, "on_editing_amplicon_pct", 10.0) - 20.0) < 1e-9
    assert 1.9 < simulated_assay_value(results, "off_editing_amplicon_pct", 10.0) < 2.1


def test_delivery_and_expression_observables_are_available():
    results = _fake_results()
    assert simulated_assay_value(results, "vector_qpcr_blood", 10.0) == 0.5
    assert simulated_assay_value(results, "puf_mrna_qpcr", 10.0) == 0.2
    assert simulated_assay_value(results, "editor_western", 10.0) == 0.3
