from __future__ import annotations

import math


def _rates_per_minute(section: dict) -> dict:
    converted: dict[str, float | str] = {}
    for key, value in section.items():
        if key.startswith("k_"):
            converted[f"{key}_per_min"] = float(value) / 60.0
        else:
            converted[key] = value
    return converted


def build_spatial_parameters(config: dict) -> dict:
    """Build the explicit minute-based parameter contract consumed by C++."""

    intracellular = _rates_per_minute(config["intracellular"])
    editing = _rates_per_minute(config["editing"])
    bbb = _rates_per_minute(config["bbb"])
    base_uptake = float(config["intracellular"]["k_ISF_to_cell"]) / 60.0

    return {
        "simulation": {
            "max_time_min": 4320.0,
            "diffusion_dt_min": 0.5,
            "mechanics_dt_min": 6.0,
            "phenotype_dt_min": 6.0,
            "intracellular_dt_min": 1.0,
            "output_interval_min": 120.0,
        },
        "microenvironment": {
            "diffusion_coefficient_um2_per_min": 10.0,
            "decay_rate_per_min": math.log(2.0) / 1440.0,
            "mesh_spacing_um": 10.0,
            "perivascular_shell_thickness_um": 10.0,
        },
        "bbb": bbb,
        "intracellular": intracellular,
        "editing": editing,
        "cell_types": {
            "endothelial": {"uptake_rate_per_min": 0.0, "apoe_scale": 0.0},
            "neuron": {"uptake_rate_per_min": base_uptake, "apoe_scale": 0.2},
            "astrocyte": {
                "uptake_rate_per_min": 1.4 * base_uptake,
                "apoe_scale": 1.0,
            },
        },
    }

