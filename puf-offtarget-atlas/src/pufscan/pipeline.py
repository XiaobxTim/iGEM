from __future__ import annotations

import hashlib
import json
import logging
import time
from collections.abc import Callable
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from pufscan import __version__
from pufscan.config import ScanConfig
from pufscan.consequence import PotentialEditingEvent, find_potential_editing_events
from pufscan.coordinates import MappedHit, TranscriptCoordinateIndex
from pufscan.expression import ExpressionMatrix, ExpressionSummary, prepare_expression
from pufscan.gencode import strip_version
from pufscan.gene_metadata import (
    EnsemblRestGeneMetadataProvider,
    GeneMetadata,
    GeneMetadataProvider,
    LocalGeneMetadataProvider,
)
from pufscan.report import generate_report
from pufscan.schemas import RunMetadata, RunResult
from pufscan.scoring import (
    SequenceScoringModel,
    SubstitutionMatrixScorer,
    UniformMismatchScorer,
    calculate_risk_score,
    load_position_weights,
    load_substitution_matrix,
)
from pufscan.search import SearchEngine, SearchHit, iter_transcripts
from pufscan.structure import RNAplfoldAnalyzer, StructureResult

LOGGER = logging.getLogger(__name__)

BINDING_SCORES = {
    "splice_junction_spanning": 1.0,
    "splice_proximal": 0.9,
    "CDS": 0.8,
    "5UTR": 0.65,
    "3UTR": 0.6,
    "exonic_non_coding": 0.5,
    "other": 0.3,
}
EDITING_SCORES = {
    "stop_gained": 1.0,
    "stop_lost": 0.95,
    "start_lost": 0.95,
    "missense": 0.85,
    "splice_proximal": 0.85,
    "splice_junction_spanning": 0.85,
    "synonymous": 0.35,
    "UTR_edit": 0.55,
    "ncRNA_edit": 0.5,
    "unknown": 0.4,
}

RESULT_COLUMNS = tuple(
    """
    rank risk_score risk_priority query_rna query_dna_normalized motif_length
    matched_sequence_rna matched_sequence_dna match_pattern exact_match mismatch_count
    mismatch_positions weighted_mismatch_penalty sequence_identity sequence_score search_orientation
    gene_id gene_id_without_version gene_name gene_type gene_description transcript_id
    transcript_id_without_version transcript_type transcript_tags chromosome strand genomic_blocks
    transcript_start transcript_end transcript_region exon_number junction_spanning
    distance_to_nearest_junction overlaps_start_codon overlaps_stop_codon target_tissue_tpm
    target_tissue_tpm_by_tissue all_tissue_tpm max_tissue_tpm median_tissue_tpm
    top_expressed_tissue number_of_tissues_tpm_ge_1 number_of_tissues_tpm_ge_10
    tissue_specificity_tau expression_score motif_mean_unpaired_probability
    motif_min_unpaired_probability opening_energy accessibility_score accessibility_profile
    structure_window_transcript_start transcript_features transcript_junctions editor editing_window
    editable_base_count potential_edit_positions coding_consequence consequence_summary
    consequence_evidence_level consequence_score missing_features warnings
    """.split()
)

ProgressCallback = Callable[[str, float], None]


def _progress(callback: ProgressCallback | None, message: str, fraction: float) -> None:
    LOGGER.info("%s (%.0f%%)", message, fraction * 100)
    if callback is not None:
        callback(message, fraction)


def _empty_expression() -> ExpressionSummary:
    return ExpressionSummary(None, None, None, None, None, None, None, None, None, None)


def _empty_structure(warning: str | None = None) -> StructureResult:
    return StructureResult(None, None, None, None, None, (), (warning,) if warning else ())


def _prepare_expression_matrix(config: ScanConfig, run_dir: Path) -> ExpressionMatrix | None:
    if config.expression is None:
        return None
    path = config.expression
    manifest_path = path.with_suffix(path.suffix + ".manifest.json")
    if path.suffix == ".parquet" and manifest_path.exists():
        return ExpressionMatrix.load(path)
    digest = hashlib.sha256(str(path.resolve()).encode()).hexdigest()[:12]
    prepared = run_dir / f"expression_{digest}.parquet"
    prepare_expression(path, prepared)
    return ExpressionMatrix.load(prepared)


def _binding_summary(region: str) -> str:
    if region in {"splice_junction_spanning", "splice_proximal"}:
        return "Potential effect on RNA processing or splicing. Experimental validation is required."
    if region == "CDS":
        return "Potential interference with translation or RNA stability. Experimental validation is required."
    return "Potential effect on RNA stability, translation, processing, or endogenous RBP competition. Experimental validation is required."


def _event_priority(event: PotentialEditingEvent) -> tuple[str, float]:
    if event.coding_consequence != "unknown":
        key = event.coding_consequence
    elif event.transcript_region in {"splice_proximal", "splice_junction_spanning"}:
        key = event.transcript_region
    elif event.transcript_region in {"5UTR", "3UTR"}:
        key = "UTR_edit"
    elif event.transcript_region == "exonic_non_coding":
        key = "ncRNA_edit"
    else:
        key = "unknown"
    return key, EDITING_SCORES[key]


def _serialize_blocks(mapped: MappedHit) -> str:
    return json.dumps(
        [{"chromosome": chrom, "start_1based": start, "end_1based": end} for chrom, start, end in mapped.genomic_blocks_1based],
        separators=(",", ":"),
    )


def _build_row(
    hit: SearchHit,
    mapped: MappedHit,
    annotation: Any,
    expression: ExpressionSummary,
    structure: StructureResult,
    events: tuple[PotentialEditingEvent, ...],
    gene_metadata: GeneMetadata,
    config: ScanConfig,
) -> dict[str, Any]:
    warnings = list(structure.warnings)
    if config.expression is None:
        warnings.append("GTEx expression data were not provided")
    if config.mode == "binding_only":
        consequence = mapped.transcript_region
        consequence_score = BINDING_SCORES.get(consequence, BINDING_SCORES["other"])
        consequence_summary = _binding_summary(consequence)
        evidence_level = "Level 1"
        coding_consequence = None
    elif events:
        prioritized = sorted(((_event_priority(event), event) for event in events), key=lambda item: item[0][1], reverse=True)
        (consequence, consequence_score), representative = prioritized[0]
        consequence_summary = representative.consequence_summary
        evidence_level = representative.consequence_evidence_level
        coding_consequence = representative.coding_consequence
    else:
        consequence = "unknown"
        consequence_score = EDITING_SCORES["unknown"]
        consequence_summary = "No potential editable base was found in the configured window."
        evidence_level = "Level 1"
        coding_consequence = "unknown"
    components = {
        "sequence": hit.sequence_score,
        "accessibility": structure.accessibility_score,
        "expression": expression.expression_score,
        "consequence": consequence_score,
    }
    risk = calculate_risk_score(components, config.risk_weights.model_dump())
    return {
        "rank": 0,
        "risk_score": risk.risk_score,
        "risk_priority": risk.risk_priority,
        "query_rna": hit.query_rna,
        "query_dna_normalized": hit.query_dna_normalized,
        "motif_length": hit.motif_length,
        "matched_sequence_rna": hit.matched_sequence_rna,
        "matched_sequence_dna": hit.matched_sequence_dna,
        "match_pattern": hit.match_pattern,
        "exact_match": hit.exact_match,
        "mismatch_count": hit.mismatch_count,
        "mismatch_positions": json.dumps(hit.mismatch_positions_1based),
        "weighted_mismatch_penalty": hit.weighted_mismatch_penalty,
        "sequence_identity": hit.sequence_identity,
        "sequence_score": hit.sequence_score,
        "search_orientation": hit.search_orientation,
        "gene_id": annotation.gene_id or hit.gene_id,
        "gene_id_without_version": strip_version(annotation.gene_id)
        or hit.gene_id_without_version,
        "gene_name": gene_metadata.gene_name or annotation.gene_name or hit.gene_name,
        "gene_type": annotation.gene_type,
        "gene_description": gene_metadata.description,
        "transcript_id": hit.transcript_id,
        "transcript_id_without_version": hit.transcript_id_without_version,
        "transcript_type": annotation.transcript_type or hit.transcript_type,
        "transcript_tags": json.dumps(annotation.tags),
        "chromosome": mapped.chromosome,
        "strand": mapped.strand,
        "genomic_blocks": _serialize_blocks(mapped),
        "transcript_start": mapped.transcript_start_1based,
        "transcript_end": mapped.transcript_end_1based,
        "transcript_region": mapped.transcript_region,
        "exon_number": json.dumps(mapped.exon_numbers),
        "junction_spanning": mapped.junction_spanning,
        "distance_to_nearest_junction": mapped.distance_to_nearest_junction,
        "overlaps_start_codon": mapped.overlaps_start_codon,
        "overlaps_stop_codon": mapped.overlaps_stop_codon,
        "target_tissue_tpm": expression.target_tissue_tpm,
        "target_tissue_tpm_by_tissue": json.dumps(expression.target_tissue_tpm_by_tissue),
        "all_tissue_tpm": json.dumps(expression.all_tissue_tpm),
        "max_tissue_tpm": expression.max_tissue_tpm,
        "median_tissue_tpm": expression.median_tissue_tpm,
        "top_expressed_tissue": expression.top_expressed_tissue,
        "number_of_tissues_tpm_ge_1": expression.number_of_tissues_tpm_ge_1,
        "number_of_tissues_tpm_ge_10": expression.number_of_tissues_tpm_ge_10,
        "tissue_specificity_tau": expression.tissue_specificity_tau,
        "expression_score": expression.expression_score,
        "motif_mean_unpaired_probability": structure.motif_mean_unpaired_probability,
        "motif_min_unpaired_probability": structure.motif_min_unpaired_probability,
        "opening_energy": structure.motif_opening_energy,
        "accessibility_score": structure.accessibility_score,
        "accessibility_profile": json.dumps(structure.local_profile),
        "structure_window_transcript_start": max(1, mapped.transcript_start_1based - config.structure.flank_nt),
        "transcript_features": json.dumps(annotation.feature_intervals),
        "transcript_junctions": json.dumps(annotation.junctions),
        "editor": config.editor,
        "editing_window": f"{config.editing_window[0]}:{config.editing_window[1]}" if config.editing_window else None,
        "editable_base_count": len(events),
        "potential_edit_positions": json.dumps([event.transcript_position_1based for event in events]),
        "coding_consequence": coding_consequence,
        "consequence_summary": consequence_summary,
        "consequence_evidence_level": evidence_level,
        "consequence_score": consequence_score,
        "missing_features": json.dumps(risk.missing_features),
        "warnings": json.dumps(warnings),
    }


def _write_bed_files(frame: pd.DataFrame, run_dir: Path) -> None:
    with (run_dir / "candidates.bed").open("w", encoding="utf-8") as bed, (run_dir / "candidates.bed12").open("w", encoding="utf-8") as bed12:
        for row in frame.to_dict(orient="records"):
            blocks = json.loads(row["genomic_blocks"])
            name = f"{row['transcript_id']}:{row['transcript_start']}-{row['transcript_end']}"
            score = max(0, min(1000, round(float(row["risk_score"]) * 10)))
            for block_index, block in enumerate(blocks, 1):
                bed.write(f"{block['chromosome']}\t{block['start_1based'] - 1}\t{block['end_1based']}\t{name}:block{block_index}\t{score}\t{row['strand']}\n")
            chrom_start0 = min(block["start_1based"] - 1 for block in blocks)
            chrom_end0 = max(block["end_1based"] for block in blocks)
            ordered = sorted(blocks, key=lambda block: block["start_1based"])
            sizes = ",".join(str(block["end_1based"] - block["start_1based"] + 1) for block in ordered) + ","
            starts = ",".join(str(block["start_1based"] - 1 - chrom_start0) for block in ordered) + ","
            bed12.write(
                f"{row['chromosome']}\t{chrom_start0}\t{chrom_end0}\t{name}\t{score}\t{row['strand']}\t{chrom_start0}\t{chrom_end0}\t8,127,140\t{len(ordered)}\t{sizes}\t{starts}\n"
            )


def _write_outputs(frame: pd.DataFrame, events: list[dict[str, Any]], run_dir: Path) -> None:
    frame.to_parquet(run_dir / "all_transcript_hits.parquet", index=False)
    frame.to_csv(run_dir / "all_transcript_hits.tsv.gz", sep="\t", index=False, compression="gzip")
    frame.head(100).to_csv(run_dir / "top_hits.tsv", sep="\t", index=False)
    event_columns = list(PotentialEditingEvent.__dataclass_fields__)
    pd.DataFrame(events, columns=event_columns).to_csv(run_dir / "potential_editing_events.tsv", sep="\t", index=False)
    if frame.empty:
        loci = pd.DataFrame(columns=["chromosome", "strand", "genomic_blocks", "transcripts", "genes", "max_risk_score", "affects_mane_select", "affects_principal"])
    else:
        loci = (
            frame.groupby(["chromosome", "strand", "genomic_blocks"], dropna=False)
            .agg(
                transcripts=("transcript_id", lambda values: json.dumps(sorted(set(values)))),
                genes=("gene_name", lambda values: json.dumps(sorted(set(str(value) for value in values)))),
                max_risk_score=("risk_score", "max"),
                affects_mane_select=("transcript_tags", lambda values: any("MANE_Select" in value for value in values)),
                affects_principal=("transcript_tags", lambda values: any("appris_principal" in value for value in values)),
            )
            .reset_index()
        )
    loci.to_csv(run_dir / "unique_genomic_loci.tsv", sep="\t", index=False)
    _write_bed_files(frame, run_dir)


def run_scan(config: ScanConfig, progress: ProgressCallback | None = None) -> RunResult:
    started = time.perf_counter()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_dir = config.output_dir / f"{config.query.replace('T', 'U')}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    metadata = RunMetadata(
        software_version=__version__,
        gencode_release=config.gencode_release,
        status="running",
        parameters=config.model_dump(mode="json"),
        input_files={"fasta": str(config.fasta), "gtf": str(config.gtf), "expression": str(config.expression) if config.expression else None},
        generated_at=datetime.now(UTC),
    )
    metadata_path = run_dir / "run_metadata.json"
    metadata_path.write_text(metadata.model_dump_json(indent=2), encoding="utf-8")
    try:
        weights = load_position_weights(config.position_weights, len(config.query)) if config.position_weights else None
        scorer: SequenceScoringModel
        if config.substitution_matrix:
            scorer = SubstitutionMatrixScorer(
                load_substitution_matrix(config.substitution_matrix, len(config.query)), weights
            )
        else:
            scorer = UniformMismatchScorer(weights)
        engine = SearchEngine(
            config.query,
            config.max_mismatches,
            config.search_reverse_complement,
            scorer,
        )
        _progress(progress, "Searching transcriptome", 0.2)
        search_result = engine.search(iter_transcripts(config.fasta))
        hit_transcripts = {hit.transcript_id for hit in search_result.hits}
        _progress(progress, "Loading candidate annotations", 0.35)
        coordinate_index = TranscriptCoordinateIndex.from_path(config.gtf, hit_transcripts)
        sequences = {record.transcript_id: record.sequence for record in iter_transcripts(config.fasta) if record.transcript_id in hit_transcripts}
        expression_matrix = _prepare_expression_matrix(config, run_dir)
        metadata_provider: GeneMetadataProvider | None = None
        if config.gene_metadata:
            metadata_provider = LocalGeneMetadataProvider(config.gene_metadata)
        elif config.ensembl_rest:
            metadata_provider = EnsemblRestGeneMetadataProvider(config.cache_dir / "ensembl")
        structure_analyzer = RNAplfoldAnalyzer(
            cache_dir=config.cache_dir / "structure",
            window_size=config.structure.window_size,
            max_base_pair_span=config.structure.max_base_pair_span,
            unpaired_length=config.structure.unpaired_length,
            temperature_c=config.structure.temperature_c,
        )
        structure_keys = {
            (hit.transcript_id, hit.transcript_start_1based, hit.search_orientation)
            for hit in sorted(search_result.hits, key=lambda item: item.sequence_score, reverse=True)[: config.structure.top_n]
        }
        rows: list[dict[str, Any]] = []
        event_rows: list[dict[str, Any]] = []
        structure_started = time.perf_counter()
        _progress(progress, "Annotating candidates", 0.45)
        for hit in search_result.hits:
            start0 = hit.transcript_start_1based - 1
            end0 = hit.transcript_end_1based
            mapped = coordinate_index.map_hit(hit.transcript_id, start0, end0, config.splice_proximity_nt)
            annotation = coordinate_index.get(hit.transcript_id)
            stable_gene = str(strip_version(annotation.gene_id) or hit.gene_id_without_version or "")
            gene_metadata = (
                metadata_provider.get_gene_metadata(stable_gene)
                if metadata_provider
                else GeneMetadata(stable_gene, annotation.gene_name)
            )
            expression = expression_matrix.summarize(stable_gene, config.target_tissues) if expression_matrix else _empty_expression()
            sequence = sequences[hit.transcript_id]
            key = (hit.transcript_id, hit.transcript_start_1based, hit.search_orientation)
            if config.structure.enabled and key in structure_keys:
                local_start = max(0, start0 - config.structure.flank_nt)
                local_end = min(len(sequence), end0 + config.structure.flank_nt)
                structure = structure_analyzer.analyze(sequence[local_start:local_end], start0 - local_start, end0 - local_start)
            elif config.structure.enabled:
                structure = _empty_structure("Structure not calculated because candidate was outside structure_top_n")
            else:
                structure = _empty_structure()
            events: tuple[PotentialEditingEvent, ...] = ()
            if config.mode == "editor_fusion" and config.editor and config.editing_window:
                events = find_potential_editing_events(
                    hit.transcript_id,
                    sequence,
                    start0,
                    end0,
                    config.editor,
                    config.editing_window,
                    coordinate_index,
                    config.splice_proximity_nt,
                )
                event_rows.extend(asdict(event) for event in events)
            rows.append(
                _build_row(
                    hit,
                    mapped,
                    annotation,
                    expression,
                    structure,
                    events,
                    gene_metadata,
                    config,
                )
            )
        structure_time = time.perf_counter() - structure_started
        frame = pd.DataFrame(rows, columns=RESULT_COLUMNS)
        if not frame.empty:
            frame = frame.sort_values(["risk_score", "sequence_score"], ascending=False).reset_index(drop=True)
            frame["rank"] = range(1, len(frame) + 1)
        _progress(progress, "Writing outputs", 0.78)
        _write_outputs(frame, event_rows, run_dir)
        summary = {
            "input_motif": config.query.replace("T", "U"),
            "motif_length": len(config.query),
            "maximum_mismatches": config.max_mismatches,
            "gencode_release": config.gencode_release,
            "candidate_sites": len(frame),
            "candidate_genes": int(frame["gene_id"].nunique()) if not frame.empty else 0,
            "candidate_transcripts": int(frame["transcript_id"].nunique()) if not frame.empty else 0,
            "unique_genomic_loci": int(frame[["chromosome", "strand", "genomic_blocks"]].drop_duplicates().shape[0]) if not frame.empty else 0,
            "runtime_seconds": round(time.perf_counter() - started, 6),
        }
        (run_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        generate_report(run_dir)
        runtime = time.perf_counter() - started
        metadata.status = "complete"
        metadata.runtime_seconds = runtime
        metadata.benchmark = {
            "transcripts_scanned": search_result.statistics.transcripts_scanned,
            "total_nucleotides_scanned": search_result.statistics.total_nucleotides_scanned,
            "candidate_count": search_result.statistics.candidate_count,
            "search_wall_time_seconds": search_result.statistics.search_wall_time_seconds,
            "structure_wall_time_seconds": structure_time,
            "peak_memory_mb": None,
        }
        metadata_path.write_text(metadata.model_dump_json(indent=2), encoding="utf-8")
        _progress(progress, "Complete", 1.0)
        return RunResult(output_dir=run_dir, candidate_count=len(frame), summary=summary)
    except Exception as error:
        metadata.status = "failed"
        metadata.runtime_seconds = time.perf_counter() - started
        metadata.error = str(error)
        metadata_path.write_text(metadata.model_dump_json(indent=2), encoding="utf-8")
        raise
