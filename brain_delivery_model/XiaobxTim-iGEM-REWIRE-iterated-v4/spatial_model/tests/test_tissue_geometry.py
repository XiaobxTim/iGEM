from __future__ import annotations

import csv

import numpy as np
from scipy.spatial.distance import pdist, squareform

from spatial_model.tissue_geometry import (
    TissueGeometry,
    generate_tissue,
    write_cells_csv,
)


def test_default_tissue_has_exact_cell_counts_and_vessel_layout():
    geometry = TissueGeometry()

    cells = generate_tissue(geometry, seed=42)

    names, counts = np.unique(cells.cell_type, return_counts=True)
    assert dict(zip(names, counts, strict=True)) == {
        "astrocyte": 400,
        "endothelial": 280,
        "neuron": 800,
    }
    endothelial = cells.select("endothelial")
    vessel_radius = np.hypot(endothelial.y, endothelial.z)
    assert np.allclose(vessel_radius, geometry.vessel_radius_um)
    assert endothelial.x.min() >= geometry.x_min_um
    assert endothelial.x.max() <= geometry.x_max_um


def test_default_tissue_is_reproducible_nonoverlapping_and_inside_domain():
    geometry = TissueGeometry()
    first = generate_tissue(geometry, seed=42)
    second = generate_tissue(geometry, seed=42)

    assert np.array_equal(first.xyz, second.xyz)
    assert np.array_equal(first.cell_type, second.cell_type)
    assert np.all(first.x - first.radius_um >= geometry.x_min_um)
    assert np.all(first.x + first.radius_um <= geometry.x_max_um)
    assert np.all(first.y - first.radius_um >= geometry.y_min_um)
    assert np.all(first.y + first.radius_um <= geometry.y_max_um)
    assert np.all(first.z - first.radius_um >= geometry.z_min_um)
    assert np.all(first.z + first.radius_um <= geometry.z_max_um)

    distances = squareform(pdist(first.xyz))
    required = first.radius_um[:, None] + first.radius_um[None, :]
    np.fill_diagonal(distances, np.inf)
    assert np.min(distances - required) >= -1e-9

    parenchyma = first.cell_type != "endothelial"
    radial = np.hypot(first.y[parenchyma], first.z[parenchyma])
    required_clearance = (
        geometry.vessel_radius_um
        + geometry.endothelial_radius_um
        + first.radius_um[parenchyma]
    )
    assert np.all(radial >= required_clearance)


def test_cells_csv_has_physicell_and_paraview_columns(tmp_path):
    cells = generate_tissue(TissueGeometry(neuron_count=3, astrocyte_count=2), seed=7)
    output = tmp_path / "cells.csv"

    write_cells_csv(cells, output)

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "cell_id",
        "cell_type",
        "x",
        "y",
        "z",
        "radius_um",
        "distance_to_vessel_um",
    ]
    assert len(rows) == 285

