from __future__ import annotations

import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path


def parse_lunp(path: Path) -> dict[int, tuple[float | None, ...]]:
    table: dict[int, tuple[float | None, ...]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            position = int(fields[0])
            values = tuple(None if value.upper() == "NA" else float(value) for value in fields[1:])
            table[position] = values
    return table


@dataclass(frozen=True)
class StructureResult:
    motif_mean_unpaired_probability: float | None
    motif_min_unpaired_probability: float | None
    motif_unpaired_probability: float | None
    motif_opening_energy: float | None
    accessibility_score: float | None
    local_profile: tuple[float | None, ...]
    warnings: tuple[str, ...] = ()


class RNAplfoldAnalyzer:
    def __init__(
        self,
        cache_dir: Path,
        window_size: int = 150,
        max_base_pair_span: int = 100,
        unpaired_length: int = 12,
        temperature_c: float = 37.0,
        executable: str = "RNAplfold",
    ):
        self.cache_dir = cache_dir
        self.window_size = window_size
        self.max_base_pair_span = max_base_pair_span
        self.unpaired_length = unpaired_length
        self.temperature_c = temperature_c
        self.executable = executable

    def analyze(self, sequence: str, motif_start0: int, motif_end0: int) -> StructureResult:
        motif_length = motif_end0 - motif_start0
        if motif_length <= 0 or motif_start0 < 0 or motif_end0 > len(sequence):
            raise ValueError("Motif interval lies outside local sequence")
        if motif_length > self.unpaired_length:
            raise ValueError("RNAplfold unpaired_length must cover the motif length")
        executable = shutil.which(self.executable)
        if executable is None:
            return StructureResult(None, None, None, None, None, (), ("RNAplfold is not available; structure fields are NA",))
        cache_payload = {
            "sequence": sequence.upper().replace("T", "U"),
            "motif_start0": motif_start0,
            "motif_end0": motif_end0,
            "window_size": self.window_size,
            "max_base_pair_span": self.max_base_pair_span,
            "unpaired_length": self.unpaired_length,
            "temperature_c": self.temperature_c,
        }
        digest = hashlib.sha256(json.dumps(cache_payload, sort_keys=True).encode()).hexdigest()
        cache_path = self.cache_dir / f"{digest}.json"
        if cache_path.exists():
            data = json.loads(cache_path.read_text(encoding="utf-8"))
            data["local_profile"] = tuple(data["local_profile"])
            data["warnings"] = tuple(data["warnings"])
            return StructureResult(**data)
        with tempfile.TemporaryDirectory(prefix="pufscan-plfold-") as temporary:
            completed = subprocess.run(
                [
                    executable,
                    "-W",
                    str(min(self.window_size, len(sequence))),
                    "-L",
                    str(min(self.max_base_pair_span, len(sequence))),
                    "-u",
                    str(self.unpaired_length),
                    "-T",
                    str(self.temperature_c),
                ],
                input=str(cache_payload["sequence"]) + "\n",
                text=True,
                cwd=temporary,
                capture_output=True,
                check=False,
            )
            if completed.returncode != 0:
                warning = f"RNAplfold failed: {completed.stderr.strip() or 'unknown error'}"
                return StructureResult(None, None, None, None, None, (), (warning,))
            candidates = list(Path(temporary).glob("*_lunp"))
            if not candidates:
                return StructureResult(None, None, None, None, None, (), ("RNAplfold produced no _lunp file",))
            table = parse_lunp(candidates[0])
        profile = tuple(table.get(position, (None,))[0] for position in range(1, len(sequence) + 1))
        motif_values = [value for value in profile[motif_start0:motif_end0] if value is not None]
        whole_values = table.get(motif_start0 + 1, ())
        whole_probability = whole_values[motif_length - 1] if len(whole_values) >= motif_length else None
        opening_energy = None
        if whole_probability is not None and whole_probability > 0:
            gas_constant_kcal = 0.0019872041
            opening_energy = -gas_constant_kcal * (273.15 + self.temperature_c) * math.log(whole_probability)
        result = StructureResult(
            motif_mean_unpaired_probability=sum(motif_values) / len(motif_values) if motif_values else None,
            motif_min_unpaired_probability=min(motif_values) if motif_values else None,
            motif_unpaired_probability=whole_probability,
            motif_opening_energy=opening_energy,
            accessibility_score=sum(motif_values) / len(motif_values) if motif_values else None,
            local_profile=profile,
        )
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
        return result
