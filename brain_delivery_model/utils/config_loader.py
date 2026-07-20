from __future__ import annotations

from pathlib import Path
from typing import Dict, Any
import yaml


def load_yaml(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data or {}


def load_base_config(project_root: str | Path) -> Dict[str, Any]:
    root = Path(project_root)
    return load_yaml(root / "config" / "base_config.yaml")
