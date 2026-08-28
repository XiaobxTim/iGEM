from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import loadmat


OUTPUT_PATTERN = re.compile(r"output\d+\.xml$")


@dataclass(frozen=True)
class SnapshotDescriptor:
    xml_path: Path
    time_min: float


@dataclass(frozen=True)
class Snapshot:
    time_min: float
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    field: np.ndarray
    positions: np.ndarray
    cell_data: dict[str, np.ndarray]


def _required(root: ET.Element, path: str) -> ET.Element:
    node = root.find(path)
    if node is None:
        raise ValueError(f"MultiCellDS snapshot is missing required element: {path}")
    return node


def _float_vector(node: ET.Element) -> np.ndarray:
    return np.fromstring(node.text or "", sep=" ", dtype=float)


def discover_snapshots(raw_dir: str | Path) -> list[SnapshotDescriptor]:
    """Return numbered PhysiCell snapshots in chronological order."""

    directory = Path(raw_dir)
    descriptors: list[SnapshotDescriptor] = []
    for xml_path in sorted(directory.glob("output*.xml")):
        if not OUTPUT_PATTERN.fullmatch(xml_path.name):
            continue
        root = ET.parse(xml_path).getroot()
        current_time = _required(root, "metadata/current_time")
        descriptors.append(
            SnapshotDescriptor(xml_path=xml_path, time_min=float(current_time.text or "nan"))
        )

    if not descriptors:
        raise FileNotFoundError(f"no numbered PhysiCell snapshots found in {directory}")
    descriptors.sort(key=lambda item: (item.time_min, item.xml_path.name))
    times = [item.time_min for item in descriptors]
    if not np.all(np.isfinite(times)) or len(times) != len(set(times)):
        raise ValueError(f"snapshot times must be finite and unique: {times}")
    return descriptors


def _load_mat_matrix(path: Path, key: str) -> np.ndarray:
    contents = loadmat(path)
    if key not in contents:
        raise ValueError(f"{path} does not contain MATLAB variable {key!r}")
    matrix = np.asarray(contents[key], dtype=float)
    if matrix.ndim != 2 or not np.all(np.isfinite(matrix)):
        raise ValueError(f"{path}:{key} must be a finite 2-D matrix")
    return matrix


def load_snapshot(descriptor: SnapshotDescriptor) -> Snapshot:
    """Load a PhysiCell MultiCellDS XML/MAT snapshot into NumPy arrays."""

    root = ET.parse(descriptor.xml_path).getroot()
    base = descriptor.xml_path.parent
    mesh = _required(root, "microenvironment/domain/mesh")
    x = _float_vector(_required(mesh, "x_coordinates"))
    y = _float_vector(_required(mesh, "y_coordinates"))
    z = _float_vector(_required(mesh, "z_coordinates"))
    if min(x.size, y.size, z.size) == 0:
        raise ValueError(f"empty mesh coordinates in {descriptor.xml_path}")

    environment_filename = (
        _required(root, "microenvironment/domain/data/filename").text or ""
    ).strip()
    environment = _load_mat_matrix(
        base / environment_filename, "multiscale_microenvironment"
    )
    voxel_count = x.size * y.size * z.size
    if environment.shape[1] != voxel_count or environment.shape[0] < 5:
        raise ValueError(
            f"unexpected microenvironment shape {environment.shape}; "
            f"expected at least (5, {voxel_count})"
        )
    # MultiCellDS stores x, y, z, voxel volume, then one row per density.
    field = environment[4].reshape((z.size, y.size, x.size), order="C")
    if np.any(field < -1e-14):
        raise ValueError("extracellular_AAV contains negative concentrations")
    field = np.maximum(field, 0.0)

    simplified = _required(
        root,
        "cellular_information/cell_populations/cell_population/custom/simplified_data",
    )
    cell_filename = (_required(simplified, "filename").text or "").strip()
    cells = _load_mat_matrix(base / cell_filename, "cells")
    labels = _required(simplified, "labels")
    cell_data: dict[str, np.ndarray] = {}
    positions: np.ndarray | None = None
    for label in labels.findall("label"):
        name = (label.text or "").strip()
        index = int(label.attrib["index"])
        size = int(label.attrib.get("size", "1"))
        if index < 0 or size < 1 or index + size > cells.shape[0]:
            raise ValueError(f"invalid cell label {name!r} at row {index} size {size}")
        values = cells[index : index + size].T.copy()
        if name == "position":
            if size != 3:
                raise ValueError("PhysiCell position label must contain three rows")
            positions = values
        elif size == 1:
            # PhysiCell currently emits elapsed_time_in_phase twice; the later
            # instance represents the same scalar and safely replaces the first.
            cell_data[name] = values[:, 0]
        else:
            for component, suffix in enumerate(("x", "y", "z")[:size]):
                cell_data[f"{name}_{suffix}"] = values[:, component]

    if positions is None:
        raise ValueError(f"cell positions are missing from {descriptor.xml_path}")
    if positions.shape[0] != cells.shape[1]:
        raise ValueError("cell position count does not match cell matrix")

    return Snapshot(
        time_min=descriptor.time_min,
        x=x,
        y=y,
        z=z,
        field=field,
        positions=positions,
        cell_data=cell_data,
    )
