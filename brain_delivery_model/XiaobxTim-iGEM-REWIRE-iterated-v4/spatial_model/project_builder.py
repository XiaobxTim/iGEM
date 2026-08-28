from __future__ import annotations

import copy
import math
import xml.etree.ElementTree as ET
from pathlib import Path

from spatial_model.intracellular_reference import STATE_NAMES


CUSTOM_FIELDS = (
    *STATE_NAMES,
    "editing_fraction",
    "off_target_burden",
    "distance_to_vessel_um",
    "uptake_rate",
    "endothelial_surface_AAV",
    "endothelial_internalized_AAV",
    "BBB_release_rate",
    "vector_cell_loss_cumulative",
    "vector_nuclear_loss_cumulative",
)


def _number(value: float) -> str:
    return format(float(value), ".12g")


def _required(root: ET.Element, path: str) -> ET.Element:
    node = root.find(path)
    if node is None:
        raise ValueError(f"PhysiCell template is missing required element: {path}")
    return node


def _set_text(root: ET.Element, path: str, value: object) -> None:
    _required(root, path).text = str(value)


def _set_zero_death_and_motility(cell_definition: ET.Element) -> None:
    for rate in cell_definition.findall("phenotype/death/model/death_rate"):
        rate.text = "0"
    motility = _required(cell_definition, "phenotype/motility/options/enabled")
    motility.text = "false"
    for duration in cell_definition.findall("phenotype/cycle//duration"):
        duration.text = "1e12"


def _set_cell_volume(cell_definition: ET.Element, radius_um: float) -> None:
    volume = 4.0 * math.pi * radius_um**3 / 3.0
    _required(cell_definition, "phenotype/volume/total").text = _number(volume)
    _required(cell_definition, "phenotype/volume/nuclear").text = _number(0.2 * volume)


def _set_substrate(cell_definition: ET.Element, uptake_rate: float) -> None:
    substrate = _required(cell_definition, "phenotype/secretion/substrate")
    substrate.attrib["name"] = "extracellular_AAV"
    _required(substrate, "uptake_rate").text = _number(uptake_rate)
    _required(substrate, "secretion_rate").text = "0"
    _required(substrate, "net_export_rate").text = "0"
    _required(
        cell_definition,
        "phenotype/motility/options/chemotaxis/substrate",
    ).text = "extracellular_AAV"
    for sensitivity in cell_definition.findall(
        "phenotype/motility/options/advanced_chemotaxis/chemotactic_sensitivities/chemotactic_sensitivity"
    ):
        sensitivity.attrib["substrate"] = "extracellular_AAV"


def _set_custom_data(cell_definition: ET.Element, cell_type_code: int) -> None:
    custom = _required(cell_definition, "custom_data")
    custom.clear()
    code = ET.SubElement(
        custom,
        "cell_type_code",
        {"conserved": "false", "units": "dimensionless", "description": "ParaView cell type code"},
    )
    code.text = str(cell_type_code)
    for field in CUSTOM_FIELDS:
        node = ET.SubElement(
            custom,
            field,
            {"conserved": "false", "units": "normalized", "description": "spatial brain-delivery state"},
        )
        node.text = "0"


def _build_cell_definitions(root: ET.Element, parameters: dict) -> None:
    definitions = _required(root, "cell_definitions")
    source = _required(definitions, "cell_definition")
    definitions.clear()
    specs = (
        ("default", 0, 7.5, 0.0),
        ("endothelial", 1, 5.0, 0.0),
        (
            "neuron",
            2,
            7.5,
            parameters["cell_types"]["neuron"]["uptake_rate_per_min"],
        ),
        (
            "astrocyte",
            3,
            8.5,
            parameters["cell_types"]["astrocyte"]["uptake_rate_per_min"],
        ),
    )
    for name, identifier, radius, uptake in specs:
        node = copy.deepcopy(source)
        node.attrib["name"] = name
        node.attrib["ID"] = str(identifier)
        if name == "default":
            node.attrib.pop("parent_type", None)
        else:
            node.attrib["parent_type"] = "default"
        _set_zero_death_and_motility(node)
        _set_cell_volume(node, radius)
        _set_substrate(node, float(uptake))
        _set_custom_data(node, identifier)
        definitions.append(node)


def _replace_user_parameters(
    root: ET.Element,
    boundary_csv: Path,
    parameter_file: Path,
    cells_csv: Path,
) -> None:
    user = _required(root, "user_parameters")
    user.clear()
    values = {
        "boundary_csv": ("string", str(boundary_csv)),
        "parameter_file": ("string", str(parameter_file)),
        "cells_csv": ("string", str(cells_csv)),
        "source_scale": ("double", "1"),
        "vessel_radius_um": ("double", "12"),
        "endothelial_radius_um": ("double", "5"),
        "perivascular_shell_thickness_um": ("double", "10"),
        "intracellular_dt_min": ("double", "1"),
        "random_seed": ("int", "42"),
    }
    for name, (kind, value) in values.items():
        node = ET.SubElement(
            user,
            name,
            {"type": kind, "units": "none", "description": "brain-delivery project parameter"},
        )
        node.text = value


def render_physicell_config(
    template_path: str | Path,
    output_path: str | Path,
    *,
    boundary_csv: str | Path,
    parameter_file: str | Path,
    cells_csv: str | Path,
    output_dir: str | Path,
    parameters: dict,
    max_time_min: float,
    output_interval_min: float,
    seed: int,
) -> Path:
    """Render a complete PhysiCell 1.14.2/Studio-compatible project XML."""

    template = Path(template_path)
    if not template.exists():
        raise FileNotFoundError(template)
    if max_time_min <= 0.0 or output_interval_min <= 0.0:
        raise ValueError("simulation duration and output interval must be positive")

    tree = ET.parse(template)
    root = tree.getroot()
    domain_values = {
        "x_min": -200,
        "x_max": 200,
        "y_min": -150,
        "y_max": 150,
        "z_min": -150,
        "z_max": 150,
        "dx": 10,
        "dy": 10,
        "dz": 10,
        "use_2D": "false",
    }
    for key, value in domain_values.items():
        _set_text(root, f"domain/{key}", value)
    _set_text(root, "overall/max_time", _number(max_time_min))
    _set_text(root, "overall/dt_diffusion", "0.5")
    _set_text(root, "overall/dt_mechanics", "6")
    _set_text(root, "overall/dt_phenotype", "6")
    _set_text(root, "parallel/omp_num_threads", "4")
    _set_text(root, "save/folder", str(Path(output_dir)))
    _set_text(root, "save/full_data/interval", _number(output_interval_min))
    _set_text(root, "save/full_data/enable", "true")
    _set_text(root, "save/SVG/enable", "false")
    _set_text(root, "options/random_seed", str(seed))

    variable = _required(root, "microenvironment_setup/variable")
    variable.attrib["name"] = "extracellular_AAV"
    variable.attrib["units"] = "normalized"
    _required(variable, "physical_parameter_set/diffusion_coefficient").text = _number(
        parameters["microenvironment"]["diffusion_coefficient_um2_per_min"]
    )
    _required(variable, "physical_parameter_set/decay_rate").text = _number(
        parameters["microenvironment"]["decay_rate_per_min"]
    )
    _required(variable, "initial_condition").text = "0"
    _required(variable, "Dirichlet_boundary_condition").attrib["enabled"] = "false"
    _set_text(
        root,
        "microenvironment_setup/options/track_internalized_substrates_in_each_agent",
        "true",
    )
    _build_cell_definitions(root, parameters)
    _set_text(root, "initial_conditions/cell_positions", "")
    _required(root, "initial_conditions/cell_positions").attrib["enabled"] = "false"
    _replace_user_parameters(
        root,
        Path(boundary_csv),
        Path(parameter_file),
        Path(cells_csv),
    )

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(tree, space="    ")
    tree.write(output, encoding="utf-8", xml_declaration=True)
    return output
