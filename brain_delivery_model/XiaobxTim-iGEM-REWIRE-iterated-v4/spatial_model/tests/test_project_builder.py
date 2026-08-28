from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from spatial_model.parameter_contract import build_spatial_parameters
from spatial_model.project_builder import render_physicell_config
from utils.config_loader import load_base_config


PHYSICELL_TEMPLATE = Path(
    "/private/tmp/PhysiCell-1.14.2/PhysiCell/sample_projects/template/config/PhysiCell_settings.xml"
)


@pytest.mark.skipif(not PHYSICELL_TEMPLATE.exists(), reason="PhysiCell 1.14.2 not staged")
def test_rendered_physicell_config_has_spatial_contract(tmp_path):
    output = tmp_path / "PhysiCell_settings.xml"
    params = build_spatial_parameters(load_base_config("."))

    render_physicell_config(
        PHYSICELL_TEMPLATE,
        output,
        boundary_csv=tmp_path / "brain_blood_input.csv",
        parameter_file=tmp_path / "spatial_parameters.cfg",
        cells_csv=tmp_path / "cells.csv",
        output_dir=tmp_path / "raw",
        parameters=params,
        max_time_min=60.0,
        output_interval_min=20.0,
        seed=42,
    )

    root = ET.parse(output).getroot()
    assert root.findtext("domain/x_min") == "-200"
    assert root.findtext("domain/z_max") == "150"
    assert root.findtext("domain/use_2D") == "false"
    assert root.findtext("overall/dt_diffusion") == "0.5"
    assert root.findtext("overall/max_time") == "60"
    variable = root.find("microenvironment_setup/variable")
    assert variable is not None
    assert variable.attrib["name"] == "extracellular_AAV"
    assert variable.findtext("physical_parameter_set/diffusion_coefficient") == "10"
    assert len(root.findall("cell_definitions/cell_definition")) == 4
    assert [node.attrib["name"] for node in root.findall("cell_definitions/cell_definition")] == [
        "default",
        "endothelial",
        "neuron",
        "astrocyte",
    ]
    custom_names = {
        node.tag
        for node in root.findall("cell_definitions/cell_definition[@name='neuron']/custom_data/*")
    }
    assert {"A_cell", "editor_protein", "editing_fraction", "distance_to_vessel_um"} <= custom_names
    for definition in root.findall("cell_definitions/cell_definition"):
        assert definition.findtext("phenotype/motility/options/chemotaxis/substrate") == "extracellular_AAV"
        sensitivities = definition.findall(
            "phenotype/motility/options/advanced_chemotaxis/chemotactic_sensitivities/chemotactic_sensitivity"
        )
        assert all(node.attrib["substrate"] == "extracellular_AAV" for node in sensitivities)
    assert root.findtext("user_parameters/boundary_csv") == str(tmp_path / "brain_blood_input.csv")
    assert root.findtext("save/folder") == str(tmp_path / "raw")


def test_render_config_rejects_missing_template(tmp_path):
    with pytest.raises(FileNotFoundError):
        render_physicell_config(
            tmp_path / "missing.xml",
            tmp_path / "out.xml",
            boundary_csv=tmp_path / "input.csv",
            parameter_file=tmp_path / "params.cfg",
            cells_csv=tmp_path / "cells.csv",
            output_dir=tmp_path / "raw",
            parameters={},
            max_time_min=60.0,
            output_interval_min=20.0,
            seed=42,
        )
