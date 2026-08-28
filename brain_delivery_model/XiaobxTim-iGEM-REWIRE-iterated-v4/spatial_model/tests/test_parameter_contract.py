from __future__ import annotations

import math

import pytest

from spatial_model.parameter_contract import build_spatial_parameters
from utils.config_loader import load_base_config


def test_spatial_parameters_convert_hour_rates_and_apply_cell_priors():
    config = load_base_config(".")

    params = build_spatial_parameters(config)

    assert params["intracellular"]["k_tx_per_min"] == pytest.approx(0.4 / 60.0)
    assert params["editing"]["k_cat_112_per_min"] == pytest.approx(1.0 / 60.0)
    assert params["microenvironment"]["decay_rate_per_min"] == pytest.approx(
        math.log(2.0) / 1440.0
    )
    assert params["cell_types"]["neuron"]["uptake_rate_per_min"] == pytest.approx(
        0.03 / 60.0
    )
    assert params["cell_types"]["astrocyte"]["uptake_rate_per_min"] == pytest.approx(
        1.4 * 0.03 / 60.0
    )
    assert params["cell_types"]["neuron"]["apoe_scale"] == 0.2
    assert params["cell_types"]["astrocyte"]["apoe_scale"] == 1.0
