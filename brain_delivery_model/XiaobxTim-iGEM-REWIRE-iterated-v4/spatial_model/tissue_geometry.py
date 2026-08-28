from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class TissueGeometry:
    x_min_um: float = -200.0
    x_max_um: float = 200.0
    y_min_um: float = -150.0
    y_max_um: float = 150.0
    z_min_um: float = -150.0
    z_max_um: float = 150.0
    mesh_spacing_um: float = 10.0
    vessel_radius_um: float = 12.0
    endothelial_radius_um: float = 5.0
    neuron_radius_um: float = 7.5
    astrocyte_radius_um: float = 8.5
    endothelial_rings: int = 40
    endothelial_per_ring: int = 7
    neuron_count: int = 800
    astrocyte_count: int = 400


@dataclass(frozen=True)
class CellTable:
    cell_id: np.ndarray
    cell_type: np.ndarray
    xyz: np.ndarray
    radius_um: np.ndarray
    distance_to_vessel_um: np.ndarray

    @property
    def x(self) -> np.ndarray:
        return self.xyz[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.xyz[:, 1]

    @property
    def z(self) -> np.ndarray:
        return self.xyz[:, 2]

    def select(self, cell_type: str) -> "CellTable":
        mask = self.cell_type == cell_type
        return CellTable(
            cell_id=self.cell_id[mask],
            cell_type=self.cell_type[mask],
            xyz=self.xyz[mask],
            radius_um=self.radius_um[mask],
            distance_to_vessel_um=self.distance_to_vessel_um[mask],
        )


def _endothelial_positions(geometry: TissueGeometry) -> np.ndarray:
    x_values = np.linspace(
        geometry.x_min_um + geometry.endothelial_radius_um,
        geometry.x_max_um - geometry.endothelial_radius_um,
        geometry.endothelial_rings,
    )
    angles = np.arange(geometry.endothelial_per_ring) * (
        2.0 * np.pi / geometry.endothelial_per_ring
    )
    positions = []
    for x in x_values:
        for angle in angles:
            positions.append(
                (
                    x,
                    geometry.vessel_radius_um * np.cos(angle),
                    geometry.vessel_radius_um * np.sin(angle),
                )
            )
    return np.asarray(positions, dtype=float)


def _sample_parenchymal_positions(
    geometry: TissueGeometry,
    rng: np.random.Generator,
    existing_xyz: list[np.ndarray],
    existing_radii: list[float],
    count: int,
    radius: float,
) -> list[np.ndarray]:
    accepted: list[np.ndarray] = []
    max_attempts = max(10000, count * 1000)
    vessel_outer_radius = geometry.vessel_radius_um + geometry.endothelial_radius_um

    for _ in range(max_attempts):
        if len(accepted) == count:
            return accepted
        candidate = np.array(
            [
                rng.uniform(geometry.x_min_um + radius, geometry.x_max_um - radius),
                rng.uniform(geometry.y_min_um + radius, geometry.y_max_um - radius),
                rng.uniform(geometry.z_min_um + radius, geometry.z_max_um - radius),
            ]
        )
        if np.hypot(candidate[1], candidate[2]) < vessel_outer_radius + radius:
            continue
        if existing_xyz:
            positions = np.vstack(existing_xyz)
            radii = np.asarray(existing_radii)
            distances = np.linalg.norm(positions - candidate, axis=1)
            if np.any(distances < radii + radius - 1e-9):
                continue
        existing_xyz.append(candidate)
        existing_radii.append(radius)
        accepted.append(candidate)
    raise RuntimeError(f"could not place {count} cells of radius {radius:g} um")


def generate_tissue(geometry: TissueGeometry | None = None, seed: int = 42) -> CellTable:
    """Create deterministic, non-overlapping endothelial and brain agents."""

    spec = geometry or TissueGeometry()
    if spec.endothelial_per_ring < 3 or spec.endothelial_rings < 1:
        raise ValueError("endothelial layout requires at least one ring of three cells")
    rng = np.random.default_rng(seed)

    endothelial = _endothelial_positions(spec)
    xyz_list = [position.copy() for position in endothelial]
    radii_list = [spec.endothelial_radius_um] * len(endothelial)
    types = ["endothelial"] * len(endothelial)

    neurons = _sample_parenchymal_positions(
        spec,
        rng,
        xyz_list,
        radii_list,
        spec.neuron_count,
        spec.neuron_radius_um,
    )
    types.extend(["neuron"] * len(neurons))

    astrocytes = _sample_parenchymal_positions(
        spec,
        rng,
        xyz_list,
        radii_list,
        spec.astrocyte_count,
        spec.astrocyte_radius_um,
    )
    types.extend(["astrocyte"] * len(astrocytes))

    xyz = np.vstack(xyz_list)
    radii = np.asarray(radii_list, dtype=float)
    cell_types = np.asarray(types, dtype="U16")
    vessel_outer_radius = spec.vessel_radius_um + spec.endothelial_radius_um
    distance = np.maximum(np.hypot(xyz[:, 1], xyz[:, 2]) - vessel_outer_radius, 0.0)
    return CellTable(
        cell_id=np.arange(len(xyz), dtype=np.int64),
        cell_type=cell_types,
        xyz=xyz,
        radius_um=radii,
        distance_to_vessel_um=distance,
    )


def _format_number(value: float) -> str:
    return format(float(value), ".12g")


def write_cells_csv(cells: CellTable, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "cell_id",
                "cell_type",
                "x",
                "y",
                "z",
                "radius_um",
                "distance_to_vessel_um",
            )
        )
        for index in range(len(cells.cell_id)):
            writer.writerow(
                (
                    int(cells.cell_id[index]),
                    str(cells.cell_type[index]),
                    _format_number(cells.x[index]),
                    _format_number(cells.y[index]),
                    _format_number(cells.z[index]),
                    _format_number(cells.radius_um[index]),
                    _format_number(cells.distance_to_vessel_um[index]),
                )
            )
    return output

