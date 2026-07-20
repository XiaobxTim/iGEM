from pathlib import Path

import pytest
from pydantic import ValidationError

from pufscan.config import ScanConfig, load_yaml


def test_editor_fusion_requires_editor_and_window(tmp_path: Path) -> None:
    with pytest.raises(ValidationError):
        ScanConfig(query="AACGUCUAUA", fasta=tmp_path / "x.fa", gtf=tmp_path / "x.gtf", mode="editor_fusion")


def test_binding_only_rejects_editor_parameters(tmp_path: Path) -> None:
    with pytest.raises(ValidationError):
        ScanConfig(
            query="AACGUCUAUA",
            fasta=tmp_path / "x.fa",
            gtf=tmp_path / "x.gtf",
            editor="APOBEC_C2U",
            editing_window=(-15, 10),
        )


def test_default_yaml_loads() -> None:
    config = load_yaml(Path("configs/default.yaml"))
    assert config["software"]["gencode_release"] == 50

