from __future__ import annotations

import math

import pytest

from spatial_model.parameter_contract import build_spatial_parameters, write_parameter_file
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


def test_parameter_file_is_flat_stable_and_cpp_readable(tmp_path):
    config = load_base_config(".")
    output = tmp_path / "spatial_parameters.cfg"

    write_parameter_file(build_spatial_parameters(config), output)

    lines = output.read_text(encoding="utf-8").splitlines()
    assert "intracellular.k_tx_per_min=0.00666666666667" in lines
    assert "editing.editor_type=A3A" in lines
    assert lines == sorted(lines)
