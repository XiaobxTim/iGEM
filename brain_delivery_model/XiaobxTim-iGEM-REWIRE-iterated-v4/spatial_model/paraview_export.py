from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

from spatial_model.raw_output import Snapshot, discover_snapshots, load_snapshot


def _spacing(coordinates: np.ndarray, name: str) -> float:
    if coordinates.size < 2:
        return 1.0
    differences = np.diff(coordinates)
    if np.any(differences <= 0.0) or not np.allclose(differences, differences[0]):
        raise ValueError(f"{name} coordinates must be uniformly increasing")
    return float(differences[0])


def compute_snapshot_metrics(snapshot: Snapshot) -> dict[str, float]:
    """Compute stable summary metrics used by dose comparison plots."""

    dx = _spacing(snapshot.x, "x")
    dy = _spacing(snapshot.y, "y")
    dz = _spacing(snapshot.z, "z")
    cell_types = snapshot.cell_data["cell_type"]
    editing = snapshot.cell_data["editing_fraction"]
    off_target = snapshot.cell_data["off_target_burden"]
    neurons = cell_types == 2
    astrocytes = cell_types == 3
    brain_cells = neurons | astrocytes

    def mean(values: np.ndarray, mask: np.ndarray) -> float:
        return float(np.mean(values[mask])) if np.any(mask) else 0.0

    release = snapshot.cell_data.get("BBB_release_rate", np.zeros_like(cell_types))
    endothelial = cell_types == 1
    return {
        "time_min": float(snapshot.time_min),
        "extracellular_AAV_mass": float(snapshot.field.sum() * dx * dy * dz),
        "neuron_editing_mean": mean(editing, neurons),
        "astrocyte_editing_mean": mean(editing, astrocytes),
        "brain_cell_editing_mean": mean(editing, brain_cells),
        "brain_cell_off_target_mean": mean(off_target, brain_cells),
        "BBB_release_rate": mean(release, endothelial),
    }


def write_pvd(
    path: str | Path, entries: list[tuple[float, Path | str]]
) -> Path:
    """Write a ParaView collection manifest with deterministic ordering."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    root = ET.Element(
        "VTKFile",
        {"type": "Collection", "version": "0.1", "byte_order": "LittleEndian"},
    )
    collection = ET.SubElement(root, "Collection")
    times: set[float] = set()
    for time_min, filename in sorted(entries, key=lambda item: item[0]):
        time_value = float(time_min)
        if not np.isfinite(time_value) or time_value in times:
            raise ValueError("PVD timesteps must be finite and unique")
        times.add(time_value)
        ET.SubElement(
            collection,
            "DataSet",
            {
                "timestep": format(time_value, ".12g"),
                "group": "",
                "part": "0",
                "file": Path(filename).as_posix(),
            },
        )
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    return output


def align_scenario_series(
    scenarios: list[tuple[str, str | Path]],
) -> list[tuple[float, list[tuple[str, Path]]]]:
    """Resolve several PVD files and require a common comparison time grid."""

    resolved: list[tuple[str, list[tuple[float, Path]]]] = []
    for label, filename in scenarios:
        pvd = Path(filename).resolve()
        datasets = ET.parse(pvd).getroot().findall("Collection/DataSet")
        entries = [
            (float(node.attrib["timestep"]), pvd.parent / node.attrib["file"])
            for node in datasets
        ]
        if not entries:
            raise ValueError(f"empty scenario series: {pvd}")
        resolved.append((label, entries))
    if not resolved:
        raise ValueError("at least one scenario is required")
    reference = [time for time, _ in resolved[0][1]]
    if any([time for time, _ in entries] != reference for _, entries in resolved[1:]):
        raise ValueError("comparison scenarios must have identical time grids")
    return [
        (
            time,
            [(label, entries[index][1]) for label, entries in resolved],
        )
        for index, time in enumerate(reference)
    ]


def write_comparison_series(
    scenarios: list[tuple[str, str | Path]], output_dir: str | Path
) -> Path:
    """Create a ParaView multiblock series with scenarios translated side-by-side."""

    try:
        import vtk
    except ImportError as error:
        raise RuntimeError("comparison export must be run with ParaView pvpython") from error

    aligned = align_scenario_series(scenarios)
    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    frames = output / "frames"
    frames.mkdir()
    offsets = (np.arange(len(scenarios), dtype=float) - (len(scenarios) - 1) / 2.0) * 450.0
    pvd_entries: list[tuple[float, Path]] = []
    for frame_index, (time_min, inputs) in enumerate(aligned):
        comparison = vtk.vtkMultiBlockDataSet()
        comparison.SetNumberOfBlocks(len(inputs))
        for scenario_index, ((label, vtm_path), x_offset) in enumerate(zip(inputs, offsets, strict=True)):
            reader = vtk.vtkXMLMultiBlockDataReader()
            reader.SetFileName(str(vtm_path))
            reader.Update()
            source = reader.GetOutput()
            if source.GetNumberOfBlocks() != 2:
                raise ValueError(f"scenario frame must contain field and cell blocks: {vtm_path}")

            shifted = vtk.vtkMultiBlockDataSet()
            shifted.SetNumberOfBlocks(2)
            image = vtk.vtkImageData()
            image.ShallowCopy(source.GetBlock(0))
            origin = image.GetOrigin()
            image.SetOrigin(origin[0] + float(x_offset), origin[1], origin[2])
            shifted.SetBlock(0, image)
            shifted.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(), "microenvironment")

            transform = vtk.vtkTransform()
            transform.Translate(float(x_offset), 0.0, 0.0)
            cell_filter = vtk.vtkTransformPolyDataFilter()
            cell_filter.SetTransform(transform)
            cell_filter.SetInputData(source.GetBlock(1))
            cell_filter.Update()
            cells = vtk.vtkPolyData()
            cells.ShallowCopy(cell_filter.GetOutput())
            shifted.SetBlock(1, cells)
            shifted.GetMetaData(1).Set(vtk.vtkCompositeDataSet.NAME(), "cells")
            comparison.SetBlock(scenario_index, shifted)
            comparison.GetMetaData(scenario_index).Set(vtk.vtkCompositeDataSet.NAME(), label)

        frame_name = Path("frames") / f"comparison_{frame_index:04d}.vtm"
        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetFileName(str(output / frame_name))
        writer.SetInputData(comparison)
        writer.SetCompressorTypeToZLib()
        if writer.Write() != 1:
            raise RuntimeError(f"failed to write comparison frame {frame_name}")
        pvd_entries.append((time_min, frame_name))
    return write_pvd(output / "comparison.pvd", pvd_entries)


def _write_vtm(path: Path, field_file: Path, cell_file: Path) -> Path:
    root = ET.Element(
        "VTKFile",
        {"type": "vtkMultiBlockDataSet", "version": "1.0", "byte_order": "LittleEndian"},
    )
    multiblock = ET.SubElement(root, "vtkMultiBlockDataSet")
    ET.SubElement(
        multiblock,
        "DataSet",
        {"index": "0", "name": "microenvironment", "file": field_file.as_posix()},
    )
    ET.SubElement(
        multiblock,
        "DataSet",
        {"index": "1", "name": "cells", "file": cell_file.as_posix()},
    )
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
    return path


def _add_point_array(dataset, numpy_support, values: np.ndarray, name: str) -> None:
    contiguous = np.ascontiguousarray(values)
    vtk_array = numpy_support.numpy_to_vtk(contiguous, deep=True)
    vtk_array.SetName(name)
    dataset.GetPointData().AddArray(vtk_array)


def _write_field_vti(snapshot: Snapshot, path: Path) -> Path:
    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError as error:
        raise RuntimeError("VTK export must be run with ParaView pvpython") from error

    nx, ny, nz = snapshot.x.size, snapshot.y.size, snapshot.z.size
    dx = _spacing(snapshot.x, "x")
    dy = _spacing(snapshot.y, "y")
    dz = _spacing(snapshot.z, "z")
    image = vtk.vtkImageData()
    image.SetDimensions(nx, ny, nz)
    image.SetOrigin(float(snapshot.x[0]), float(snapshot.y[0]), float(snapshot.z[0]))
    image.SetSpacing(dx, dy, dz)

    concentration = snapshot.field.ravel(order="C")
    zz, yy, xx = np.meshgrid(snapshot.z, snapshot.y, snapshot.x, indexing="ij")
    radial = np.hypot(yy, zz)
    distance = np.maximum(radial - 17.0, 0.0).ravel(order="C")
    vessel_mask = (radial < 17.0).astype(np.uint8).ravel(order="C")
    tissue_mask = (radial >= 17.0).astype(np.uint8).ravel(order="C")
    shell = ((radial >= 17.0) & (radial < 27.0)).ravel(order="C")
    cell_types = snapshot.cell_data["cell_type"]
    endothelial = cell_types == 1
    release = snapshot.cell_data.get("BBB_release_rate", np.zeros_like(cell_types))
    release_value = float(np.mean(release[endothelial])) if np.any(endothelial) else 0.0
    release_field = np.where(shell, release_value, 0.0)

    _add_point_array(image, numpy_support, concentration, "extracellular_AAV")
    _add_point_array(
        image,
        numpy_support,
        np.log10(np.maximum(concentration, 1e-18)),
        "log10_extracellular_AAV",
    )
    _add_point_array(image, numpy_support, distance, "distance_to_vessel_um")
    _add_point_array(image, numpy_support, vessel_mask, "vessel_mask")
    _add_point_array(image, numpy_support, tissue_mask, "tissue_mask")
    _add_point_array(image, numpy_support, release_field, "BBB_release_rate")
    image.GetPointData().SetActiveScalars("extracellular_AAV")

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(image)
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToAppended()
    if writer.Write() != 1:
        raise RuntimeError(f"failed to write {path}")
    return path


def _write_cells_vtp(snapshot: Snapshot, path: Path) -> Path:
    try:
        import vtk
        from vtk.util import numpy_support
    except ImportError as error:
        raise RuntimeError("VTK export must be run with ParaView pvpython") from error

    points = vtk.vtkPoints()
    points.SetData(numpy_support.numpy_to_vtk(np.ascontiguousarray(snapshot.positions), deep=True))
    vertices = vtk.vtkCellArray()
    for index in range(snapshot.positions.shape[0]):
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(index)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetVerts(vertices)
    source_types = snapshot.cell_data["cell_type"].astype(np.int32)
    # Compact rendering code: endothelial=0, neuron=1, astrocyte=2.
    mapped_types = np.where(source_types > 0, source_types - 1, -1).astype(np.int32)
    _add_point_array(polydata, numpy_support, mapped_types, "cell_type")
    _add_point_array(polydata, numpy_support, source_types, "physicell_cell_type")
    for name, values in sorted(snapshot.cell_data.items()):
        export_name = {"ID": "cell_id", "radius": "radius_um"}.get(name, name)
        if export_name in {"cell_type", "physicell_cell_type"}:
            continue
        array = np.asarray(values)
        if array.ndim == 1 and array.shape[0] == snapshot.positions.shape[0]:
            if export_name == "cell_id":
                array = array.astype(np.int64)
            _add_point_array(polydata, numpy_support, array, export_name)
    polydata.GetPointData().SetActiveScalars("cell_type")

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetCompressorTypeToZLib()
    writer.SetDataModeToAppended()
    if writer.Write() != 1:
        raise RuntimeError(f"failed to write {path}")
    return path


def publish_atomically(staging: str | Path, destination: str | Path) -> Path:
    """Rename a complete staged directory into place without overwriting data."""

    source = Path(staging)
    target = Path(destination)
    if not source.is_dir():
        raise FileNotFoundError(source)
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, target)
    return target


def export_scenario(raw_dir: str | Path, output_dir: str | Path) -> Path:
    """Convert all raw snapshots and atomically publish one ParaView series."""

    raw = Path(raw_dir).resolve()
    target = Path(output_dir).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}-", dir=target.parent))
    try:
        fields_dir = staging / "fields"
        cells_dir = staging / "cells"
        frames_dir = staging / "frames"
        for directory in (fields_dir, cells_dir, frames_dir):
            directory.mkdir()

        entries: list[tuple[float, Path]] = []
        metrics: list[dict[str, float]] = []
        for frame_index, descriptor in enumerate(discover_snapshots(raw)):
            snapshot = load_snapshot(descriptor)
            stem = f"frame_{frame_index:04d}"
            field_name = Path("fields") / f"{stem}.vti"
            cell_name = Path("cells") / f"{stem}.vtp"
            frame_name = Path("frames") / f"{stem}.vtm"
            _write_field_vti(snapshot, staging / field_name)
            _write_cells_vtp(snapshot, staging / cell_name)
            _write_vtm(
                staging / frame_name,
                Path("..") / field_name,
                Path("..") / cell_name,
            )
            entries.append((snapshot.time_min, frame_name))
            metrics.append(compute_snapshot_metrics(snapshot))

        write_pvd(staging / "simulation.pvd", entries)
        with (staging / "metrics.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
            writer.writeheader()
            writer.writerows(metrics)
        (staging / "README.txt").write_text(
            "Open simulation.pvd in ParaView 6.0.1. The microenvironment block "
            "contains volume fields and the cells block contains point arrays.\n",
            encoding="utf-8",
        )
        return publish_atomically(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert PhysiCell output to ParaView VTK")
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = _parse_args()
    result = export_scenario(arguments.raw_dir, arguments.output_dir)
    print(f"ParaView series written to {result / 'simulation.pvd'}")


if __name__ == "__main__":
    main()
