from __future__ import annotations

from pathlib import Path

from spatial_model.physicell_stage import stage_physicell_project


PROJECT_ROOT = Path(__file__).resolve().parents[2]
STAGED_PHYSICELL = Path("/private/tmp/PhysiCell-1.14.2/PhysiCell")


def test_stage_physicell_project_installs_overlay_and_build_contract():
    staged = stage_physicell_project(STAGED_PHYSICELL, PROJECT_ROOT / "spatial_model")

    assert staged == STAGED_PHYSICELL
    assert (staged / "main.cpp").is_file()
    assert (staged / "custom_modules" / "brain_delivery.cpp").is_file()
    assert (staged / "custom_modules" / "intracellular_model.cpp").is_file()
    makefile = (staged / "Makefile").read_text(encoding="utf-8")
    assert "PROGRAM_NAME := brain_delivery" in makefile
    assert "brain_delivery.o intracellular_model.o" in makefile
    assert "-Xpreprocessor -fopenmp" in makefile
