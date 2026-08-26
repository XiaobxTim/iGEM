from __future__ import annotations

from models.full_model.simulator import build_initial_state
from models.full_model.state_vector import STATE_ORDER, build_index_map
from utils.config_loader import load_base_config


def test_config_supports_all_administration_routes() -> None:
    config = load_base_config(".")

    assert {"footpad", "im", "iv"} <= set(config["absorption"])


def test_iv_dose_starts_in_blood_instead_of_local_depot() -> None:
    config = load_base_config(".")
    config["route"] = "iv"
    idx = build_index_map(STATE_ORDER)

    initial = build_initial_state(2.5, config)

    assert initial[idx["A_dep"]] == 0.0
    assert initial[idx["A_blood"]] == 2.5
