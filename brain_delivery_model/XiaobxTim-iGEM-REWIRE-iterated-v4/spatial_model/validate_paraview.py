from __future__ import annotations

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


FIELD_ARRAYS = {
    "extracellular_AAV",
    "log10_extracellular_AAV",
    "distance_to_vessel_um",
    "vessel_mask",
    "tissue_mask",
    "BBB_release_rate",
}
CELL_ARRAYS = {
    "cell_id",
    "cell_type",
    "physicell_cell_type",
    "radius_um",
    "editor_protein",
    "editing_fraction",
    "off_target_burden",
}


def _array_names(attributes) -> set[str]:
    return {
        attributes.GetArrayName(index)
        for index in range(attributes.GetNumberOfArrays())
    }


def validate_series(pvd_path: str | Path) -> dict[str, object]:
    """Round-trip a generated series through the VTK shipped with ParaView."""

    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError as error:
        raise RuntimeError("validation must be run with ParaView pvpython") from error

    pvd = Path(pvd_path).resolve()
    root = ET.parse(pvd).getroot()
    datasets = root.findall("Collection/DataSet")
    if not datasets:
        raise ValueError(f"PVD collection is empty: {pvd}")
    times: list[float] = []
    cell_counts: list[int] = []
    dimensions: list[tuple[int, int, int]] = []
    for dataset in datasets:
        times.append(float(dataset.attrib["timestep"]))
        vtm = pvd.parent / dataset.attrib["file"]
        vtm_root = ET.parse(vtm).getroot()
        blocks = {
            node.attrib["name"]: vtm.parent / node.attrib["file"]
            for node in vtm_root.findall("vtkMultiBlockDataSet/DataSet")
        }
        if set(blocks) != {"microenvironment", "cells"}:
            raise ValueError(f"unexpected VTM blocks in {vtm}: {sorted(blocks)}")

        image_reader = vtk.vtkXMLImageDataReader()
        image_reader.SetFileName(str(blocks["microenvironment"]))
        image_reader.Update()
        image = image_reader.GetOutput()
        field_names = _array_names(image.GetPointData())
        if not FIELD_ARRAYS <= field_names:
            raise ValueError(f"missing field arrays: {sorted(FIELD_ARRAYS - field_names)}")
        concentration = numpy_support.vtk_to_numpy(
            image.GetPointData().GetArray("extracellular_AAV")
        )
        if not np.all(np.isfinite(concentration)) or np.any(concentration < 0.0):
            raise ValueError("extracellular_AAV must be finite and non-negative")
        dimensions.append(tuple(image.GetDimensions()))

        cell_reader = vtk.vtkXMLPolyDataReader()
        cell_reader.SetFileName(str(blocks["cells"]))
        cell_reader.Update()
        cells = cell_reader.GetOutput()
        cell_names = _array_names(cells.GetPointData())
        if not CELL_ARRAYS <= cell_names:
            raise ValueError(f"missing cell arrays: {sorted(CELL_ARRAYS - cell_names)}")
        if cells.GetNumberOfVerts() != cells.GetNumberOfPoints():
            raise ValueError("each exported cell must have one vertex")
        cell_counts.append(cells.GetNumberOfPoints())

    if times != sorted(times) or len(times) != len(set(times)):
        raise ValueError("PVD times must be sorted and unique")
    return {
        "pvd": str(pvd),
        "frames": len(datasets),
        "time_min": times,
        "dimensions": sorted(set(dimensions)),
        "cell_counts": sorted(set(cell_counts)),
        "status": "valid",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate ParaView brain-delivery data")
    parser.add_argument("pvd", type=Path)
    arguments = parser.parse_args()
    print(json.dumps(validate_series(arguments.pvd), indent=2))


if __name__ == "__main__":
    main()
