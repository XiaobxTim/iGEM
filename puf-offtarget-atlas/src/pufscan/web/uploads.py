from __future__ import annotations

import shutil
import uuid
from pathlib import Path

from fastapi import UploadFile

from pufscan.gencode import parse_gtf
from pufscan.index import prepare_transcriptome
from pufscan.search import iter_transcripts
from pufscan.transcriptomes import TranscriptomeSpec

ALLOWED_FASTA_SUFFIXES = (".fa", ".fasta", ".fa.gz", ".fasta.gz")
ALLOWED_GTF_SUFFIXES = (".gtf", ".gtf.gz")


def _matching_suffix(filename: str | None, allowed: tuple[str, ...]) -> str:
    lowered = (filename or "").lower()
    for suffix in sorted(allowed, key=len, reverse=True):
        if lowered.endswith(suffix):
            return suffix
    raise ValueError(f"Unsupported file type: {filename or 'unnamed upload'}")


async def save_upload(upload: UploadFile, destination: Path, max_bytes: int) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    try:
        with destination.open("wb") as handle:
            while chunk := await upload.read(1024 * 1024):
                written += len(chunk)
                if written > max_bytes:
                    raise OverflowError(f"Upload exceeds the configured {max_bytes} byte limit")
                handle.write(chunk)
    except Exception:
        destination.unlink(missing_ok=True)
        raise
    finally:
        await upload.close()
    return destination


def transcript_id_overlap(fasta_path: Path, gtf_path: Path, minimum: float = 0.9) -> float:
    fasta_ids = {record.transcript_id for record in iter_transcripts(fasta_path)}
    annotation_ids = {
        record.transcript_id for record in parse_gtf(gtf_path) if record.transcript_id is not None
    }
    if not fasta_ids:
        raise ValueError("The uploaded FASTA contains no transcript records")
    overlap = len(fasta_ids & annotation_ids) / len(fasta_ids)
    if overlap < minimum:
        raise ValueError(
            f"FASTA and GTF transcript identifiers overlap by {overlap:.1%}; at least 90% is required"
        )
    return overlap


def prepare_custom_upload(
    staging_dir: Path,
    final_dir: Path,
    identifier: str,
    display_name: str,
    species: str,
    assembly: str,
    provider: str,
    release: str,
    fasta_path: Path,
    gtf_path: Path,
) -> TranscriptomeSpec:
    transcript_id_overlap(fasta_path, gtf_path)
    prepared_dir = staging_dir / "prepared"
    prepare_transcriptome(
        fasta_path,
        gtf_path,
        prepared_dir,
        species=species,
        assembly=assembly,
        provider=provider,
        release=release,
    )
    if final_dir.exists():
        raise ValueError(f"Transcriptome storage for {identifier!r} already exists")
    staging_dir.replace(final_dir)
    return TranscriptomeSpec(
        id=identifier,
        display_name=display_name,
        species=species,
        assembly=assembly,
        provider=provider,
        release=release,
        fasta_path=final_dir / "prepared/transcripts.fa",
        annotation_path=final_dir / "prepared/annotation.parquet",
        source="custom",
    )


def new_staging_directory(custom_data_dir: Path) -> Path:
    path = custom_data_dir / f".staging-{uuid.uuid4().hex}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def discard_staging(path: Path) -> None:
    if path.name.startswith(".staging-") and path.exists():
        shutil.rmtree(path)


__all__ = [
    "ALLOWED_FASTA_SUFFIXES",
    "ALLOWED_GTF_SUFFIXES",
    "discard_staging",
    "new_staging_directory",
    "prepare_custom_upload",
    "save_upload",
    "transcript_id_overlap",
]
