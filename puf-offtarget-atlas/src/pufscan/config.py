from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, Field, model_validator

from pufscan.sequence import normalize_query


class StructureConfig(BaseModel):
    enabled: bool = True
    flank_nt: int = Field(default=100, ge=0)
    window_size: int = Field(default=150, ge=1)
    max_base_pair_span: int = Field(default=100, ge=1)
    unpaired_length: int = Field(default=12, ge=12)
    temperature_c: float = 37.0
    top_n: int = Field(default=5000, ge=0)


class RiskWeights(BaseModel):
    sequence: float = Field(default=0.5, ge=0)
    accessibility: float = Field(default=0.2, ge=0)
    expression: float = Field(default=0.2, ge=0)
    consequence: float = Field(default=0.1, ge=0)


class ScanConfig(BaseModel):
    query: str
    fasta: Path
    gtf: Path
    expression: Path | None = None
    target_tissues: tuple[str, ...] = ()
    max_mismatches: int = Field(default=2, ge=0, le=3)
    search_reverse_complement: bool = False
    mode: Literal["binding_only", "editor_fusion"] = "binding_only"
    editor: Literal["APOBEC_C2U", "ADAR_A2I"] | None = None
    editing_window: tuple[int, int] | None = None
    splice_proximity_nt: int = Field(default=20, ge=0)
    structure: StructureConfig = Field(default_factory=StructureConfig)
    risk_weights: RiskWeights = Field(default_factory=RiskWeights)
    output_dir: Path = Path("results")
    cache_dir: Path = Path(".pufscan_cache")
    threads: int = Field(default=1, ge=1)
    gencode_release: int = Field(default=50, ge=1)
    position_weights: Path | None = None
    substitution_matrix: Path | None = None
    gene_metadata: Path | None = None
    ensembl_rest: bool = False

    @model_validator(mode="after")
    def validate_mode(self) -> ScanConfig:
        self.query = normalize_query(self.query)
        if self.mode == "editor_fusion" and (self.editor is None or self.editing_window is None):
            raise ValueError("editor_fusion mode requires both editor and editing_window")
        if self.mode == "binding_only" and (self.editor is not None or self.editing_window is not None):
            raise ValueError("binding_only mode does not accept editor parameters")
        if self.editing_window is not None and self.editing_window[0] > self.editing_window[1]:
            raise ValueError("editing_window start must not exceed end")
        return self


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML configuration root must be a mapping")
    return data

