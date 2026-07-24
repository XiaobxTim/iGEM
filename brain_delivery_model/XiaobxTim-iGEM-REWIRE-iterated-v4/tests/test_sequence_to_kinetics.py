from __future__ import annotations

from models.editing.sequence_to_kinetics import design_to_kinetic_modifiers


def test_sequence_features_increase_intended_kinetics():
    weak = {
        "design_id": "weak",
        "editor_type": "A3A",
        "puf_repeats": 8,
        "mismatch_count": 2,
        "puf_score_112": 0.2,
        "puf_score_158": 0.1,
        "accessibility_112": 0.2,
        "accessibility_158": 0.2,
        "distance_112": 8,
        "distance_158": 8,
        "uc_context_112": 0.0,
        "uc_context_158": 0.0,
        "local_bystander_risk": 0.1,
        "deaminase_background_prior": 1.0,
    }
    strong = {
        **weak,
        "design_id": "strong",
        "puf_repeats": 10,
        "mismatch_count": 0,
        "puf_score_112": 0.9,
        "accessibility_112": 0.8,
        "distance_112": 2,
        "uc_context_112": 1.0,
    }

    weak_mod = design_to_kinetic_modifiers(weak)
    strong_mod = design_to_kinetic_modifiers(strong)

    assert strong_mod["kon112_scale"] > weak_mod["kon112_scale"]
    assert strong_mod["kcat112_scale"] > weak_mod["kcat112_scale"]
    assert strong_mod["puf_offtarget_scale"] < weak_mod["puf_offtarget_scale"]
