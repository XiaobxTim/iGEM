from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import json
import numpy as np


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_results(results: Dict[str, Any], output_dir: str | Path) -> None:
    out = ensure_dir(output_dir)

    np.savez(
        out / "simulation_results.npz",
        t=results["t"],
        y=results["y"],
        state_order=np.array(results["state_order"], dtype=object),
    )

    metadata = {
        "dose": results["dose"],
        "route": results["config"].get("route", "footpad"),
        "n_states": len(results["state_order"]),
    }
    with open(out / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
