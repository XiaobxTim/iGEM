from __future__ import annotations

from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from spatial_model.paraview_export import (
    align_scenario_series,
    compute_snapshot_metrics,
    publish_atomically,
    write_pvd,
)
from spatial_model.raw_output import Snapshot


def _snapshot() -> Snapshot:
    return Snapshot(
        time_min=120.0,
        x=np.array([-5.0, 5.0]),
        y=np.array([-5.0, 5.0]),
        z=np.array([-5.0, 5.0]),
        field=np.ones((2, 2, 2)),
        positions=np.zeros((5, 3)),
        cell_data={
            "ID": np.arange(5.0),
            "cell_type": np.array([1.0, 2.0, 2.0, 3.0, 3.0]),
            "editing_fraction": np.array([0.0, 0.1, 0.3, 0.5, 0.9]),
            "off_target_burden": np.array([0.0, 0.02, 0.04, 0.06, 0.08]),
            "BBB_release_rate": np.array([0.25, 0.0, 0.0, 0.0, 0.0]),
        },
    )


def test_compute_snapshot_metrics_uses_cell_type_specific_statistics():
    metrics = compute_snapshot_metrics(_snapshot())

    assert metrics["time_min"] == 120.0
    assert metrics["extracellular_AAV_mass"] == pytest.approx(8000.0)
    assert metrics["neuron_editing_mean"] == pytest.approx(0.2)
    assert metrics["astrocyte_editing_mean"] == pytest.approx(0.7)
    assert metrics["brain_cell_off_target_mean"] == pytest.approx(0.05)
    assert metrics["BBB_release_rate"] == pytest.approx(0.25)


def test_write_pvd_records_sorted_relative_files(tmp_path):
    output = write_pvd(
        tmp_path / "simulation.pvd",
        [(10.0, Path("frames/frame_0001.vtm")), (0.0, Path("frames/frame_0000.vtm"))],
    )

    collection = ET.parse(output).getroot().find("Collection")
    datasets = collection.findall("DataSet") if collection is not None else []
    assert [float(node.attrib["timestep"]) for node in datasets] == [0.0, 10.0]
    assert [node.attrib["file"] for node in datasets] == [
        "frames/frame_0000.vtm",
        "frames/frame_0001.vtm",
    ]


def test_publish_atomically_refuses_existing_destination(tmp_path):
    staging = tmp_path / "staging"
    staging.mkdir()
    (staging / "ready.txt").write_text("ready", encoding="utf-8")
    destination = tmp_path / "published"
    destination.mkdir()

    with pytest.raises(FileExistsError):
        publish_atomically(staging, destination)


def test_align_scenario_series_rejects_different_time_grids(tmp_path):
    first = tmp_path / "first.pvd"
    second = tmp_path / "second.pvd"
    write_pvd(first, [(0.0, "a.vtm"), (10.0, "b.vtm")])
    write_pvd(second, [(0.0, "a.vtm"), (20.0, "b.vtm")])

    with pytest.raises(ValueError, match="identical time grids"):
        align_scenario_series([("low", first), ("high", second)])
