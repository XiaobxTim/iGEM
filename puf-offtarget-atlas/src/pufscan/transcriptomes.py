from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field


class TranscriptomeSpec(BaseModel):
    id: str = Field(pattern=r"^[a-z0-9][a-z0-9-]{1,62}$")
    display_name: str
    species: str
    assembly: str
    provider: str
    release: str
    fasta_path: Path
    annotation_path: Path
    expression_path: Path | None = None
    source: Literal["preset", "custom"] = "custom"


@dataclass(frozen=True)
class TranscriptomeAvailability:
    ready: bool
    missing_files: tuple[Path, ...]


def builtin_transcriptomes(project_root: Path) -> tuple[TranscriptomeSpec, ...]:
    return (
        TranscriptomeSpec(
            id="human-gencode50",
            display_name="Human · GENCODE 50",
            species="Homo sapiens",
            assembly="GRCh38.p14",
            provider="GENCODE",
            release="50",
            fasta_path=project_root / "data/gencode_v50/prepared/transcripts.fa",
            annotation_path=project_root / "data/gencode_v50/prepared/annotation.parquet",
            expression_path=project_root / "data/gtex/gene_median_tpm.parquet",
            source="preset",
        ),
        TranscriptomeSpec(
            id="mouse-gencode-m39",
            display_name="Mouse · GENCODE M39",
            species="Mus musculus",
            assembly="GRCm39",
            provider="GENCODE",
            release="M39",
            fasta_path=project_root / "data/gencode_m39/prepared/transcripts.fa",
            annotation_path=project_root / "data/gencode_m39/prepared/annotation.parquet",
            source="preset",
        ),
    )


class TranscriptomeRegistry:
    def __init__(self, database_path: Path, project_root: Path):
        self.database_path = database_path
        self.project_root = project_root
        database_path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database_path)
        connection.row_factory = sqlite3.Row
        return connection

    def _initialize(self) -> None:
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS transcriptomes (
                    id TEXT PRIMARY KEY,
                    payload TEXT NOT NULL,
                    source TEXT NOT NULL
                )
                """
            )
            for spec in builtin_transcriptomes(self.project_root):
                connection.execute(
                    "INSERT OR IGNORE INTO transcriptomes (id, payload, source) VALUES (?, ?, ?)",
                    (spec.id, spec.model_dump_json(), spec.source),
                )

    def add(self, spec: TranscriptomeSpec) -> None:
        try:
            with self._connect() as connection:
                connection.execute(
                    "INSERT INTO transcriptomes (id, payload, source) VALUES (?, ?, ?)",
                    (spec.id, spec.model_dump_json(), spec.source),
                )
        except sqlite3.IntegrityError as error:
            raise ValueError(f"Transcriptome {spec.id!r} already exists") from error

    def get(self, identifier: str) -> TranscriptomeSpec:
        with self._connect() as connection:
            row = connection.execute(
                "SELECT payload FROM transcriptomes WHERE id = ?", (identifier,)
            ).fetchone()
        if row is None:
            raise KeyError(f"Unknown transcriptome {identifier!r}")
        return TranscriptomeSpec.model_validate(json.loads(row["payload"]))

    def list(self) -> tuple[TranscriptomeSpec, ...]:
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT payload FROM transcriptomes ORDER BY source DESC, id"
            ).fetchall()
        return tuple(TranscriptomeSpec.model_validate(json.loads(row["payload"])) for row in rows)

    def availability(self, identifier: str) -> TranscriptomeAvailability:
        spec = self.get(identifier)
        required = [spec.fasta_path, spec.annotation_path]
        missing = tuple(path for path in required if not path.exists())
        return TranscriptomeAvailability(ready=not missing, missing_files=missing)


__all__ = [
    "TranscriptomeAvailability",
    "TranscriptomeRegistry",
    "TranscriptomeSpec",
    "builtin_transcriptomes",
]
