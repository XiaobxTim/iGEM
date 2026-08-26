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


def test_scan_config_records_transcriptome_metadata(tmp_path: Path) -> None:
    config = ScanConfig(
        query="AACGUCUAUA",
        fasta=tmp_path / "mouse.fa",
        gtf=tmp_path / "mouse.gtf",
        species="Mus musculus",
        genome_build="GRCm39",
        annotation_provider="GENCODE",
        annotation_release="M39",
    )

    assert config.species == "Mus musculus"
    assert config.genome_build == "GRCm39"
    assert config.annotation_release == "M39"
