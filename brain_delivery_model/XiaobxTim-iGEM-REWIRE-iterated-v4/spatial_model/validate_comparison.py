from __future__ import annotations

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def validate_comparison(pvd_path: str | Path, expected_scenarios: int = 3) -> dict[str, object]:
    try:
        import vtk
    except ImportError as error:
        raise RuntimeError("validation must be run with ParaView pvpython") from error

    pvd = Path(pvd_path).resolve()
    datasets = ET.parse(pvd).getroot().findall("Collection/DataSet")
    if not datasets:
        raise ValueError(f"empty comparison PVD: {pvd}")
    centers_by_frame: list[list[float]] = []
    for dataset in datasets:
        frame = pvd.parent / dataset.attrib["file"]
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader.SetFileName(str(frame))
        reader.Update()
        comparison = reader.GetOutput()
        if comparison.GetNumberOfBlocks() != expected_scenarios:
            raise ValueError(f"expected {expected_scenarios} scenarios in {frame}")
        centers: list[float] = []
        for index in range(expected_scenarios):
            scenario = comparison.GetBlock(index)
            if scenario is None or scenario.GetNumberOfBlocks() != 2:
                raise ValueError(f"invalid scenario block {index} in {frame}")
            image, cells = scenario.GetBlock(0), scenario.GetBlock(1)
            if image.GetNumberOfPoints() != 40 * 30 * 30 or cells.GetNumberOfPoints() != 1480:
                raise ValueError(f"unexpected data size in scenario block {index}")
            bounds = image.GetBounds()
            centers.append(0.5 * (bounds[0] + bounds[1]))
        if centers != sorted(centers) or any(
            right - left < 400.0 for left, right in zip(centers, centers[1:])
        ):
            raise ValueError(f"comparison scenarios are not separated side-by-side: {centers}")
        centers_by_frame.append(centers)
    return {
        "pvd": str(pvd),
        "frames": len(datasets),
        "scenarios": expected_scenarios,
        "x_centers_um": centers_by_frame[0],
        "status": "valid",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate side-by-side ParaView series")
    parser.add_argument("pvd", type=Path)
    arguments = parser.parse_args()
    print(json.dumps(validate_comparison(arguments.pvd), indent=2))


if __name__ == "__main__":
    main()
