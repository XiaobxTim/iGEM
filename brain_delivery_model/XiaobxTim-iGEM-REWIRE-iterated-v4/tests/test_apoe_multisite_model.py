from __future__ import annotations

import copy

from utils.config_loader import load_base_config
from models.full_model.simulator import run_simulation
from models.editing.module5 import compute_editing_metrics


def test_apoe_multisite_model_produces_site_specific_fractions():
    config = load_base_config(".")
    config = copy.deepcopy(config)
    config["editing"]["model"] = "apoe_multisite"

    results = run_simulation(config=config, dose=1.0, t_end=24.0, dt=1.0)
    metrics = compute_editing_metrics(results["y"], results["idx"], config)

    assert results["success"]
    assert "apoe4_fraction" in metrics
    assert "apoe3_like_fraction" in metrics
    assert "apoe2_like_fraction" in metrics
    assert "local_bystander_burden" in metrics
    assert metrics["apoe3_like_fraction"][-1] >= 0.0
    assert metrics["apoe2_like_fraction"][-1] >= 0.0
