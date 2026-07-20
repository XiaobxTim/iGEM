from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq

from pufscan.gencode import GtfRecord, parse_gtf, strip_version


@dataclass(frozen=True)
class CoordinateSegment:
    transcript_start0: int
    transcript_end0: int
    genomic_start1: int
    genomic_end1: int
    exon_number: str | None


@dataclass(frozen=True)
class MappedHit:
    chromosome: str
    strand: str
    transcript_start_1based: int
    transcript_end_1based: int
    genomic_blocks_1based: tuple[tuple[str, int, int], ...]
    exon_numbers: tuple[str, ...]
    junction_spanning: bool
    distance_to_nearest_junction: int | None
    transcript_region: str
    overlaps_start_codon: bool
    overlaps_stop_codon: bool


@dataclass
class TranscriptAnnotation:
    transcript_id: str
    gene_id: str | None
    gene_name: str | None
    gene_type: str | None
    transcript_name: str | None
    transcript_type: str | None
    chromosome: str
    strand: str
    tags: tuple[str, ...]
    exons: list[CoordinateSegment]
    feature_intervals: dict[str, list[tuple[int, int]]]

    @property
    def junctions(self) -> tuple[int, ...]:
        return tuple(segment.transcript_end0 for segment in self.exons[:-1])


class TranscriptCoordinateIndex:
    def __init__(self, transcripts: dict[str, TranscriptAnnotation]):
        self.transcripts = transcripts
        self._stable_ids = {strip_version(key): value for key, value in transcripts.items()}

    @classmethod
    def from_gtf(
        cls, path: Path, transcript_ids: set[str] | None = None
    ) -> TranscriptCoordinateIndex:
        records: Iterable[GtfRecord] = parse_gtf(path)
        if transcript_ids is not None:
            stable_ids = {strip_version(identifier) for identifier in transcript_ids}
            records = (
                record
                for record in records
                if record.transcript_id in transcript_ids
                or record.transcript_id_without_version in stable_ids
            )
        return cls.from_records(records)

    @classmethod
    def from_path(
        cls, path: Path, transcript_ids: set[str] | None = None
    ) -> TranscriptCoordinateIndex:
        if path.suffix != ".parquet":
            return cls.from_gtf(path, transcript_ids)
        filters = [("transcript_id", "in", sorted(transcript_ids))] if transcript_ids else None
        rows = pq.read_table(path, filters=filters).to_pylist()
        return cls.from_records(
            GtfRecord(**{**row, "tags": tuple(row.get("tags") or ())}) for row in rows
        )

    @classmethod
    def from_records(cls, records: Iterable[GtfRecord]) -> TranscriptCoordinateIndex:
        grouped: dict[str, list[GtfRecord]] = defaultdict(list)
        metadata: dict[str, GtfRecord] = {}
        for record in records:
            if record.transcript_id is None:
                continue
            grouped[record.transcript_id].append(record)
            if record.feature == "transcript":
                metadata[record.transcript_id] = record
        transcripts: dict[str, TranscriptAnnotation] = {}
        for transcript_id, records in grouped.items():
            exon_records = [record for record in records if record.feature == "exon"]
            if not exon_records:
                continue
            exemplar = metadata.get(transcript_id, exon_records[0])
            exon_records.sort(key=lambda item: item.genomic_start, reverse=exemplar.strand == "-")
            segments: list[CoordinateSegment] = []
            cursor = 0
            for exon in exon_records:
                length = exon.genomic_end - exon.genomic_start + 1
                segments.append(
                    CoordinateSegment(cursor, cursor + length, exon.genomic_start, exon.genomic_end, exon.exon_number)
                )
                cursor += length
            transcript = TranscriptAnnotation(
                transcript_id=transcript_id,
                gene_id=exemplar.gene_id,
                gene_name=exemplar.gene_name,
                gene_type=exemplar.gene_type,
                transcript_name=exemplar.transcript_name,
                transcript_type=exemplar.transcript_type,
                chromosome=exemplar.chromosome,
                strand=exemplar.strand,
                tags=exemplar.tags,
                exons=segments,
                feature_intervals=defaultdict(list),
            )
            for record in records:
                if record.feature in {"CDS", "UTR", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"}:
                    transcript.feature_intervals[record.feature].extend(
                        cls._genomic_feature_to_transcript(record, segments, exemplar.strand)
                    )
            transcripts[transcript_id] = transcript
        return cls(transcripts)

    @staticmethod
    def _genomic_feature_to_transcript(
        feature: GtfRecord, segments: list[CoordinateSegment], strand: str
    ) -> list[tuple[int, int]]:
        intervals: list[tuple[int, int]] = []
        for segment in segments:
            overlap_start = max(feature.genomic_start, segment.genomic_start1)
            overlap_end = min(feature.genomic_end, segment.genomic_end1)
            if overlap_start > overlap_end:
                continue
            if strand == "+":
                start0 = segment.transcript_start0 + overlap_start - segment.genomic_start1
                end0 = segment.transcript_start0 + overlap_end - segment.genomic_start1 + 1
            else:
                start0 = segment.transcript_start0 + segment.genomic_end1 - overlap_end
                end0 = segment.transcript_start0 + segment.genomic_end1 - overlap_start + 1
            intervals.append((start0, end0))
        return intervals

    def get(self, transcript_id: str) -> TranscriptAnnotation:
        if transcript_id in self.transcripts:
            return self.transcripts[transcript_id]
        stable = strip_version(transcript_id)
        if stable in self._stable_ids:
            return self._stable_ids[stable]
        raise KeyError(f"Transcript {transcript_id} is absent from the GTF index")

    def map_hit(
        self, transcript_id: str, start0: int, end0: int, splice_proximity_nt: int = 20
    ) -> MappedHit:
        transcript = self.get(transcript_id)
        transcript_length = transcript.exons[-1].transcript_end0
        if start0 < 0 or end0 <= start0 or end0 > transcript_length:
            raise ValueError("Hit interval lies outside transcript bounds")
        genomic_blocks: list[tuple[str, int, int]] = []
        exon_numbers: list[str] = []
        for exon in transcript.exons:
            overlap_start = max(start0, exon.transcript_start0)
            overlap_end = min(end0, exon.transcript_end0)
            if overlap_start >= overlap_end:
                continue
            if transcript.strand == "+":
                genomic_start = exon.genomic_start1 + overlap_start - exon.transcript_start0
                genomic_end = exon.genomic_start1 + overlap_end - exon.transcript_start0 - 1
            else:
                genomic_start = exon.genomic_end1 - (overlap_end - exon.transcript_start0) + 1
                genomic_end = exon.genomic_end1 - (overlap_start - exon.transcript_start0)
            genomic_blocks.append((transcript.chromosome, genomic_start, genomic_end))
            if exon.exon_number is not None:
                exon_numbers.append(exon.exon_number)
        genomic_blocks.sort(key=lambda block: block[1])
        junction_spanning = any(start0 < junction < end0 for junction in transcript.junctions)
        distance = min((min(abs(start0 - j), abs(end0 - j)) for j in transcript.junctions), default=None)
        region = self._classify_region(transcript, start0, end0, junction_spanning, distance, splice_proximity_nt)
        return MappedHit(
            chromosome=transcript.chromosome,
            strand=transcript.strand,
            transcript_start_1based=start0 + 1,
            transcript_end_1based=end0,
            genomic_blocks_1based=tuple(genomic_blocks),
            exon_numbers=tuple(exon_numbers),
            junction_spanning=junction_spanning,
            distance_to_nearest_junction=0 if junction_spanning else distance,
            transcript_region=region,
            overlaps_start_codon=self._overlaps(transcript.feature_intervals.get("start_codon", []), start0, end0),
            overlaps_stop_codon=self._overlaps(transcript.feature_intervals.get("stop_codon", []), start0, end0),
        )

    @staticmethod
    def _overlaps(intervals: list[tuple[int, int]], start0: int, end0: int) -> bool:
        return any(start0 < interval_end and end0 > interval_start for interval_start, interval_end in intervals)

    @classmethod
    def _classify_region(
        cls,
        transcript: TranscriptAnnotation,
        start0: int,
        end0: int,
        junction_spanning: bool,
        distance: int | None,
        proximity: int,
    ) -> str:
        if junction_spanning:
            return "splice_junction_spanning"
        if distance is not None and distance <= proximity:
            return "splice_proximal"
        if cls._overlaps(transcript.feature_intervals.get("CDS", []), start0, end0):
            return "CDS"
        if cls._overlaps(transcript.feature_intervals.get("five_prime_utr", []), start0, end0):
            return "5UTR"
        if cls._overlaps(transcript.feature_intervals.get("three_prime_utr", []), start0, end0):
            return "3UTR"
        if cls._overlaps(transcript.feature_intervals.get("UTR", []), start0, end0):
            cds = transcript.feature_intervals.get("CDS", [])
            if cds:
                cds_start = min(interval[0] for interval in cds)
                cds_end = max(interval[1] for interval in cds)
                if end0 <= cds_start:
                    return "5UTR"
                if start0 >= cds_end:
                    return "3UTR"
        if not transcript.feature_intervals.get("CDS"):
            return "exonic_non_coding"
        return "other"
