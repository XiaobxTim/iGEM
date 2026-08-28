from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from spatial_model.raw_output import discover_snapshots, load_snapshot


SMOKE_RAW = Path("/private/tmp/physicell-smoke.fS3jqd/raw")


@pytest.mark.skipif(not SMOKE_RAW.exists(), reason="PhysiCell smoke output not available")
def test_discover_snapshots_uses_unique_output_timepoints():
    snapshots = discover_snapshots(SMOKE_RAW)

    assert [snapshot.time_min for snapshot in snapshots] == [0.0, 5.0, 10.0]
    assert [snapshot.xml_path.name for snapshot in snapshots] == [
        "output00000000.xml",
        "output00000001.xml",
        "output00000002.xml",
    ]


@pytest.mark.skipif(not SMOKE_RAW.exists(), reason="PhysiCell smoke output not available")
def test_load_snapshot_recovers_field_cells_and_custom_data():
    descriptor = discover_snapshots(SMOKE_RAW)[-1]

    snapshot = load_snapshot(descriptor)

    assert snapshot.field.shape == (30, 30, 40)
    assert np.all(snapshot.field >= 0.0)
    assert snapshot.positions.shape == (1480, 3)
    assert snapshot.cell_data["ID"].shape == (1480,)
    assert set(np.unique(snapshot.cell_data["cell_type"])) == {1.0, 2.0, 3.0}
    assert {"editor_protein", "editing_fraction", "off_target_burden"} <= set(
        snapshot.cell_data
    )
    assert snapshot.x[0] == -195.0
    assert snapshot.x[-1] == 195.0


def test_discover_snapshots_rejects_missing_output(tmp_path):
    with pytest.raises(FileNotFoundError):
        discover_snapshots(tmp_path)
