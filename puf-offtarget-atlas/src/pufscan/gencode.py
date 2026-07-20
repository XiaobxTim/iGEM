from __future__ import annotations

import gzip
import hashlib
import json
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import requests

LOGGER = logging.getLogger(__name__)
GENCODE_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}"


def gencode_filenames(release: int) -> dict[str, str]:
    return {
        "transcript_fasta": f"gencode.v{release}.transcripts.fa.gz",
        "reference_gtf": f"gencode.v{release}.annotation.gtf.gz",
        "all_regions_gtf": f"gencode.v{release}.chr_patch_hapl_scaff.annotation.gtf.gz",
    }


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path, timeout: float = 60.0) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    mode = "ab" if existing else "wb"
    with requests.get(url, headers=headers, stream=True, timeout=timeout) as response:
        if existing and response.status_code == 200:
            mode = "wb"
        response.raise_for_status()
        with partial.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    partial.replace(destination)
    return {"url": url, "path": str(destination), "size_bytes": destination.stat().st_size, "sha256": sha256_file(destination)}


def download_gencode(release: int, output_dir: Path, include_all_regions_gtf: bool = True) -> Path:
    names = gencode_filenames(release)
    keys = ["transcript_fasta", "reference_gtf"]
    if include_all_regions_gtf:
        keys.append("all_regions_gtf")
    manifest: dict[str, object] = {"release": release, "files": {}}
    for key in keys:
        filename = names[key]
        destination = output_dir / filename
        if destination.exists():
            record = {"url": f"{GENCODE_BASE_URL.format(release=release)}/{filename}", "path": str(destination), "size_bytes": destination.stat().st_size, "sha256": sha256_file(destination), "skipped_existing": True}
        else:
            record = download_file(f"{GENCODE_BASE_URL.format(release=release)}/{filename}", destination)
        manifest["files"][key] = record  # type: ignore[index]
    manifest_path = output_dir / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest_path


def strip_version(identifier: str | None) -> str | None:
    if identifier is None:
        return None
    return identifier.split(".", 1)[0]


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_gtf_attributes(raw: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = {}
    for key, value in re.findall(r"([^\s;]+)\s+\"([^\"]*)\"\s*;", raw):
        attributes.setdefault(key, []).append(value)
    return attributes


def _first_attribute(attributes: dict[str, list[str]], key: str) -> str | None:
    values = attributes.get(key)
    return values[0] if values else None


@dataclass(frozen=True)
class GtfRecord:
    gene_id: str | None
    gene_id_without_version: str | None
    gene_name: str | None
    gene_type: str | None
    transcript_id: str | None
    transcript_id_without_version: str | None
    transcript_name: str | None
    transcript_type: str | None
    chromosome: str
    strand: str
    feature: str
    genomic_start: int
    genomic_end: int
    exon_number: str | None
    exon_id: str | None
    phase: int | None
    tags: tuple[str, ...]


def parse_gtf(path: Path) -> Iterator[GtfRecord]:
    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF line {line_number}: expected 9 columns")
            chromosome, _, feature, start, end, _, strand, phase_raw, attributes_raw = fields
            attributes = parse_gtf_attributes(attributes_raw)

            gene_id = _first_attribute(attributes, "gene_id")
            transcript_id = _first_attribute(attributes, "transcript_id")
            phase = None if phase_raw == "." else int(phase_raw)
            yield GtfRecord(
                gene_id=gene_id,
                gene_id_without_version=strip_version(gene_id),
                gene_name=_first_attribute(attributes, "gene_name"),
                gene_type=_first_attribute(attributes, "gene_type")
                or _first_attribute(attributes, "gene_biotype"),
                transcript_id=transcript_id,
                transcript_id_without_version=strip_version(transcript_id),
                transcript_name=_first_attribute(attributes, "transcript_name"),
                transcript_type=_first_attribute(attributes, "transcript_type")
                or _first_attribute(attributes, "transcript_biotype"),
                chromosome=chromosome,
                strand=strand,
                feature=feature,
                genomic_start=int(start),
                genomic_end=int(end),
                exon_number=_first_attribute(attributes, "exon_number"),
                exon_id=_first_attribute(attributes, "exon_id"),
                phase=phase,
                tags=tuple(attributes.get("tag", [])),
            )
