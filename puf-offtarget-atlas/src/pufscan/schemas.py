from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field


class RunMetadata(BaseModel):
    software_version: str
    gencode_release: int
    species: str = "Homo sapiens"
    genome_build: str = "GRCh38"
    status: Literal["running", "complete", "failed"]
    parameters: dict[str, Any]
    input_files: dict[str, str | None]
    generated_at: datetime
    runtime_seconds: float | None = None
    benchmark: dict[str, int | float | None] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: str | None = None


class RunResult(BaseModel):
    output_dir: Path
    candidate_count: int
    summary: dict[str, Any]

