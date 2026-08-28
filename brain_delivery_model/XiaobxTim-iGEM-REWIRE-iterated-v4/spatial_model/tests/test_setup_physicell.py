from __future__ import annotations

import pytest

from spatial_model.setup_physicell import PHYSICELL_VERSION, setup_physicell


def test_setup_reuses_exact_existing_version(tmp_path):
    root = tmp_path / "install/PhysiCell"
    root.mkdir(parents=True)
    (root / "VERSION.txt").write_text(PHYSICELL_VERSION, encoding="utf-8")

    assert setup_physicell(tmp_path / "install") == root


def test_setup_refuses_incomplete_existing_destination(tmp_path):
    destination = tmp_path / "install"
    destination.mkdir()

    with pytest.raises(FileExistsError):
        setup_physicell(destination)
