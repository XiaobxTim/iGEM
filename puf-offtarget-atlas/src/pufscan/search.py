from __future__ import annotations

import gzip
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import ahocorasick
from Bio import SeqIO
from pyfaidx import Fasta

from pufscan.gencode import strip_version
from pufscan.scoring import SequenceScoringModel, UniformMismatchScorer
from pufscan.sequence import (
    dna_to_rna,
    generate_variants,
    hamming_distance,
    match_pattern,
    normalize_query,
    reverse_complement,
)


@dataclass(frozen=True)
class TranscriptRecord:
    transcript_id: str
    transcript_id_without_version: str
    gene_id: str | None
    gene_id_without_version: str | None
    gene_name: str | None
    transcript_type: str | None
    sequence: str


@dataclass(frozen=True)
class SearchHit:
    query_rna: str
    query_dna_normalized: str
    motif_length: int
    transcript_id: str
    transcript_id_without_version: str
    gene_id: str | None
    gene_id_without_version: str | None
    gene_name: str | None
    transcript_type: str | None
    transcript_start_1based: int
    transcript_end_1based: int
    matched_sequence_rna: str
    matched_sequence_dna: str
    mismatch_count: int
    mismatch_positions_1based: tuple[int, ...]
    match_pattern: str
    sequence_identity: float
    sequence_score: float
    weighted_mismatch_penalty: float
    exact_match: bool
    search_orientation: str


@dataclass(frozen=True)
class SearchStatistics:
    transcripts_scanned: int
    total_nucleotides_scanned: int
    candidate_count: int
    search_wall_time_seconds: float
    variants_indexed: int


@dataclass(frozen=True)
class SearchResult:
    hits: tuple[SearchHit, ...]
    statistics: SearchStatistics


def _open_fasta(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def iter_transcripts(path: Path) -> Iterator[TranscriptRecord]:
    if path.suffix != ".gz" and Path(str(path) + ".fai").exists():
        indexed = Fasta(str(path), key_function=lambda key: key.split("|", 1)[0])
        try:
            # pyfaidx iteration yields FastaRecord objects; keys() yields stable IDs.
            for transcript_id in indexed.keys():  # noqa: SIM118
                yield TranscriptRecord(
                    transcript_id=transcript_id,
                    transcript_id_without_version=str(strip_version(transcript_id)),
                    gene_id=None,
                    gene_id_without_version=None,
                    gene_name=None,
                    transcript_type=None,
                    sequence=str(indexed[transcript_id][:]).upper().replace("U", "T"),
                )
        finally:
            indexed.close()
        return
    with _open_fasta(path) as handle:
        for fasta_record in SeqIO.parse(handle, "fasta"):  # type: ignore[no-untyped-call]
            fields = fasta_record.id.split("|")
            transcript_id = fields[0]
            gene_id = fields[1] if len(fields) > 1 and fields[1] else None
            yield TranscriptRecord(
                transcript_id=transcript_id,
                transcript_id_without_version=str(strip_version(transcript_id)),
                gene_id=gene_id,
                gene_id_without_version=strip_version(gene_id),
                gene_name=fields[6] if len(fields) > 6 and fields[6] else None,
                transcript_type=fields[8] if len(fields) > 8 and fields[8] else None,
                sequence=str(fasta_record.seq).upper().replace("U", "T"),
            )


class SearchEngine:
    def __init__(
        self,
        query: str,
        max_mismatches: int = 2,
        search_reverse_complement: bool = False,
        scorer: SequenceScoringModel | None = None,
    ):
        self.query = normalize_query(query)
        if not 0 <= max_mismatches <= 3:
            raise ValueError("max_mismatches must be between 0 and 3 in this version")
        self.max_mismatches = max_mismatches
        self.search_reverse_complement = search_reverse_complement
        self.scorer = scorer or UniformMismatchScorer()
        patterns: dict[str, list[tuple[str, str]]] = {}
        for variant in generate_variants(self.query, max_mismatches):
            patterns.setdefault(variant, []).append(("forward", self.query))
        if search_reverse_complement:
            reverse_query = reverse_complement(self.query)
            for variant in generate_variants(reverse_query, max_mismatches):
                patterns.setdefault(variant, []).append(("reverse_complement", reverse_query))
        automaton = ahocorasick.Automaton()
        for pattern, payload in patterns.items():
            automaton.add_word(pattern, (pattern, tuple(payload)))
        automaton.make_automaton()
        self.automaton = automaton
        self.variant_count = len(patterns)

    def search(self, transcripts: Iterable[TranscriptRecord]) -> SearchResult:
        started = time.perf_counter()
        hits: list[SearchHit] = []
        seen: set[tuple[str, int, int, str]] = set()
        transcript_count = 0
        nucleotide_count = 0
        motif_length = len(self.query)
        for transcript in transcripts:
            transcript_count += 1
            nucleotide_count += len(transcript.sequence)
            for end_index, (_, payloads) in self.automaton.iter(transcript.sequence):
                start0 = end_index - motif_length + 1
                end0 = end_index + 1
                observed = transcript.sequence[start0:end0]
                for orientation, reference in payloads:
                    mismatch_count, _ = hamming_distance(reference, observed)
                    if mismatch_count > self.max_mismatches:
                        continue
                    key = (transcript.transcript_id, start0, end0, orientation)
                    if key in seen:
                        continue
                    seen.add(key)
                    score = self.scorer.score(reference, observed)
                    hits.append(
                        SearchHit(
                            query_rna=dna_to_rna(self.query),
                            query_dna_normalized=self.query,
                            motif_length=motif_length,
                            transcript_id=transcript.transcript_id,
                            transcript_id_without_version=transcript.transcript_id_without_version,
                            gene_id=transcript.gene_id,
                            gene_id_without_version=transcript.gene_id_without_version,
                            gene_name=transcript.gene_name,
                            transcript_type=transcript.transcript_type,
                            transcript_start_1based=start0 + 1,
                            transcript_end_1based=end0,
                            matched_sequence_rna=dna_to_rna(observed),
                            matched_sequence_dna=observed,
                            mismatch_count=score.mismatch_count,
                            mismatch_positions_1based=score.mismatch_positions_1based,
                            match_pattern=match_pattern(reference, observed),
                            sequence_identity=score.sequence_identity,
                            sequence_score=score.sequence_score,
                            weighted_mismatch_penalty=score.weighted_mismatch_penalty,
                            exact_match=score.exact_match,
                            search_orientation=orientation,
                        )
                    )
        wall_time = time.perf_counter() - started
        hits.sort(key=lambda hit: (hit.transcript_id, hit.transcript_start_1based, hit.search_orientation))
        return SearchResult(
            hits=tuple(hits),
            statistics=SearchStatistics(
                transcripts_scanned=transcript_count,
                total_nucleotides_scanned=nucleotide_count,
                candidate_count=len(hits),
                search_wall_time_seconds=wall_time,
                variants_indexed=self.variant_count,
            ),
        )
