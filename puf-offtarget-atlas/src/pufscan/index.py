from __future__ import annotations

import gzip
import json
import shutil
from dataclasses import asdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from pyfaidx import Fasta

from pufscan.coordinates import TranscriptCoordinateIndex
from pufscan.gencode import parse_gtf, sha256_file


def _copy_or_decompress(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as incoming, destination.open("wb") as outgoing:
            shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)
    else:
        shutil.copyfile(source, destination)


def prepare_gencode(fasta_path: Path, gtf_path: Path, output_dir: Path, batch_size: int = 50000) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared_fasta = output_dir / "transcripts.fa"
    _copy_or_decompress(fasta_path, prepared_fasta)
    fasta = Fasta(str(prepared_fasta), rebuild=True, key_function=lambda key: key.split("|", 1)[0])
    transcript_count = len(fasta.keys())
    fasta.close()
    annotation_path = output_dir / "annotation.parquet"
    writer: pq.ParquetWriter | None = None
    batch: list[dict[str, object]] = []
    record_count = 0
    for record in parse_gtf(gtf_path):
        row = asdict(record)
        row["tags"] = list(record.tags)
        batch.append(row)
        if len(batch) >= batch_size:
            table = pa.Table.from_pylist(batch)
            writer = writer or pq.ParquetWriter(annotation_path, table.schema, compression="zstd")
            writer.write_table(table)
            record_count += len(batch)
            batch.clear()
    if batch:
        table = pa.Table.from_pylist(batch)
        writer = writer or pq.ParquetWriter(annotation_path, table.schema, compression="zstd")
        writer.write_table(table)
        record_count += len(batch)
    if writer is not None:
        writer.close()
    index = TranscriptCoordinateIndex.from_gtf(gtf_path)
    segment_rows: list[dict[str, object]] = []
    for transcript in index.transcripts.values():
        for exon in transcript.exons:
            segment_rows.append(
                {
                    "transcript_id": transcript.transcript_id,
                    "chromosome": transcript.chromosome,
                    "strand": transcript.strand,
                    **asdict(exon),
                }
            )
    pq.write_table(pa.Table.from_pylist(segment_rows), output_dir / "transcript_segments.parquet", compression="zstd")
    manifest = {
        "format_version": 1,
        "source_fasta": str(fasta_path),
        "source_gtf": str(gtf_path),
        "source_fasta_sha256": sha256_file(fasta_path),
        "source_gtf_sha256": sha256_file(gtf_path),
        "prepared_fasta": str(prepared_fasta),
        "fasta_index": str(prepared_fasta.with_suffix(".fa.fai")),
        "annotation_parquet": str(annotation_path),
        "transcript_segments_parquet": str(output_dir / "transcript_segments.parquet"),
        "transcript_count": transcript_count,
        "annotation_record_count": record_count,
    }
    manifest_path = output_dir / "prepared_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest_path

