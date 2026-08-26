from __future__ import annotations

import json
import sqlite3
import uuid
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

RunStatus = Literal[
    "queued",
    "preparing",
    "searching",
    "annotating",
    "scoring",
    "reporting",
    "completed",
    "failed",
]


@dataclass(frozen=True)
class RunRecord:
    id: str
    status: RunStatus
    stage: str
    progress: float
    config: dict[str, object]
    transcriptome_id: str
    output_dir: Path | None
    error: str | None
    created_at: datetime
    updated_at: datetime


class JobStore:
    def __init__(self, database_path: Path):
        self.database_path = database_path
        database_path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS runs (
                    id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    stage TEXT NOT NULL,
                    progress REAL NOT NULL,
                    config_json TEXT NOT NULL,
                    transcriptome_id TEXT NOT NULL,
                    output_dir TEXT,
                    error TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.database_path, timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    @staticmethod
    def _from_row(row: sqlite3.Row) -> RunRecord:
        return RunRecord(
            id=str(row["id"]),
            status=row["status"],
            stage=str(row["stage"]),
            progress=float(row["progress"]),
            config=json.loads(row["config_json"]),
            transcriptome_id=str(row["transcriptome_id"]),
            output_dir=Path(row["output_dir"]) if row["output_dir"] else None,
            error=row["error"],
            created_at=datetime.fromisoformat(row["created_at"]),
            updated_at=datetime.fromisoformat(row["updated_at"]),
        )

    def create(self, config: dict[str, object], transcriptome_id: str) -> RunRecord:
        identifier = uuid.uuid4().hex
        now = datetime.now(UTC).isoformat()
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO runs (
                    id, status, stage, progress, config_json, transcriptome_id,
                    output_dir, error, created_at, updated_at
                ) VALUES (?, 'queued', 'Queued', 0, ?, ?, NULL, NULL, ?, ?)
                """,
                (identifier, json.dumps(config), transcriptome_id, now, now),
            )
        return self.get(identifier)

    def get(self, identifier: str) -> RunRecord:
        with self._connect() as connection:
            row = connection.execute("SELECT * FROM runs WHERE id = ?", (identifier,)).fetchone()
        if row is None:
            raise KeyError(f"Unknown run {identifier!r}")
        return self._from_row(row)

    def claim_next(self) -> RunRecord | None:
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                "SELECT id FROM runs WHERE status = 'queued' ORDER BY created_at LIMIT 1"
            ).fetchone()
            if row is None:
                return None
            identifier = str(row["id"])
            now = datetime.now(UTC).isoformat()
            connection.execute(
                "UPDATE runs SET status = 'preparing', stage = 'Preparing analysis', updated_at = ? WHERE id = ?",
                (now, identifier),
            )
            claimed = connection.execute("SELECT * FROM runs WHERE id = ?", (identifier,)).fetchone()
        assert claimed is not None
        return self._from_row(claimed)

    def update_progress(
        self, identifier: str, stage: str, progress: float, status: RunStatus
    ) -> None:
        now = datetime.now(UTC).isoformat()
        with self._connect() as connection:
            connection.execute(
                "UPDATE runs SET stage = ?, progress = ?, status = ?, updated_at = ? WHERE id = ?",
                (stage, min(max(progress, 0.0), 1.0), status, now, identifier),
            )

    def complete(self, identifier: str, output_dir: Path) -> None:
        now = datetime.now(UTC).isoformat()
        with self._connect() as connection:
            connection.execute(
                """
                UPDATE runs SET status = 'completed', stage = 'Complete', progress = 1,
                output_dir = ?, error = NULL, updated_at = ? WHERE id = ?
                """,
                (str(output_dir), now, identifier),
            )

    def fail(self, identifier: str, error: str) -> None:
        now = datetime.now(UTC).isoformat()
        with self._connect() as connection:
            connection.execute(
                "UPDATE runs SET status = 'failed', stage = 'Analysis failed', error = ?, updated_at = ? WHERE id = ?",
                (error, now, identifier),
            )


__all__ = ["JobStore", "RunRecord", "RunStatus"]
