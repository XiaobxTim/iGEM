from pathlib import Path

import pytest

from pufscan.transcriptomes import TranscriptomeRegistry, TranscriptomeSpec


def make_spec(tmp_path: Path, identifier: str = "custom-demo") -> TranscriptomeSpec:
    fasta = tmp_path / "transcripts.fa"
    annotation = tmp_path / "annotation.parquet"
    fasta.write_text(">TX1\nAACGTCTATA\n", encoding="utf-8")
    annotation.write_bytes(b"parquet-placeholder")
    return TranscriptomeSpec(
        id=identifier,
        display_name="Demo transcriptome",
        species="Testus organismus",
        assembly="Demo1",
        provider="Custom",
        release="1",
        fasta_path=fasta,
        annotation_path=annotation,
        source="custom",
    )


def test_registry_seeds_human_and_mouse(tmp_path: Path) -> None:
    registry = TranscriptomeRegistry(tmp_path / "atlas.sqlite3", tmp_path)

    human = registry.get("human-gencode50")
    mouse = registry.get("mouse-gencode-m39")

    assert human.species == "Homo sapiens"
    assert human.assembly == "GRCh38.p14"
    assert mouse.species == "Mus musculus"
    assert mouse.assembly == "GRCm39"


def test_registry_persists_custom_transcriptome(tmp_path: Path) -> None:
    database = tmp_path / "atlas.sqlite3"
    registry = TranscriptomeRegistry(database, tmp_path)
    spec = make_spec(tmp_path)

    registry.add(spec)

    reopened = TranscriptomeRegistry(database, tmp_path)
    assert reopened.get(spec.id) == spec
    assert spec.id in {item.id for item in reopened.list()}


def test_registry_rejects_duplicate_identifier(tmp_path: Path) -> None:
    registry = TranscriptomeRegistry(tmp_path / "atlas.sqlite3", tmp_path)
    spec = make_spec(tmp_path)
    registry.add(spec)

    with pytest.raises(ValueError, match="already exists"):
        registry.add(spec)


def test_registry_reports_missing_dataset_files(tmp_path: Path) -> None:
    registry = TranscriptomeRegistry(tmp_path / "atlas.sqlite3", tmp_path)
    spec = make_spec(tmp_path)
    registry.add(spec)
    spec.fasta_path.unlink()

    availability = registry.availability(spec.id)

    assert availability.ready is False
    assert availability.missing_files == (spec.fasta_path,)
