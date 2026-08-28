import json
from pathlib import Path

import pandas as pd

from pufscan.config import ScanConfig, StructureConfig
from pufscan.index import prepare_gencode, prepare_transcriptome
from pufscan.pipeline import run_scan

DATA = Path(__file__).parent / "data"


def test_binding_only_end_to_end_writes_all_core_outputs(tmp_path: Path) -> None:
    config = ScanConfig(
        query="AACGUCUAUA",
        fasta=DATA / "synthetic.fa",
        gtf=DATA / "synthetic.gtf",
        expression=DATA / "expression.tsv",
        target_tissues=("Liver",),
        max_mismatches=1,
        structure=StructureConfig(enabled=False),
        output_dir=tmp_path,
    )
    result = run_scan(config)
    expected = {
        "run_metadata.json",
        "all_transcript_hits.parquet",
        "all_transcript_hits.tsv.gz",
        "unique_genomic_loci.tsv",
        "potential_editing_events.tsv",
        "top_hits.tsv",
        "candidates.bed",
        "candidates.bed12",
        "report.html",
        "summary.json",
        "brain_candidate_panel.csv",
        "brain_candidate_panel.metadata.json",
    }
    assert expected <= {path.name for path in result.output_dir.iterdir()}
    hits = pd.read_parquet(result.output_dir / "all_transcript_hits.parquet")
    assert len(hits) == 6
    assert {"gene_name", "transcript_region", "genomic_blocks", "risk_score"} <= set(hits.columns)
    assert hits["risk_score"].between(0, 100).all()
    brain_panel = pd.read_csv(result.output_dir / "brain_candidate_panel.csv")
    assert len(brain_panel) == len(
        hits[["chromosome", "strand", "genomic_blocks"]].drop_duplicates()
    )
    assert brain_panel["binding_score"].between(0, 1).all()
    summary = json.loads((result.output_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["candidate_sites"] == 6
    assert "prioritization score, not a calibrated probability" in (result.output_dir / "report.html").read_text()


def test_editor_fusion_generates_potential_events(tmp_path: Path) -> None:
    config = ScanConfig(
        query="AACGUCUAUA",
        fasta=DATA / "synthetic.fa",
        gtf=DATA / "synthetic.gtf",
        max_mismatches=0,
        mode="editor_fusion",
        editor="APOBEC_C2U",
        editing_window=(-2, 12),
        structure=StructureConfig(enabled=False),
        output_dir=tmp_path,
    )
    result = run_scan(config)
    events = pd.read_csv(result.output_dir / "potential_editing_events.tsv", sep="\t")
    assert len(events) > 0
    assert set(events["editor"]) == {"APOBEC_C2U"}
    assert set(events["event_label"]) == {"potential editable base"}


def test_prepared_fasta_and_annotation_index_are_reusable(tmp_path: Path) -> None:
    prepared = tmp_path / "prepared"
    prepare_gencode(DATA / "synthetic.fa", DATA / "synthetic.gtf", prepared)
    result = run_scan(
        ScanConfig(
            query="AACGUCUAUA",
            fasta=prepared / "transcripts.fa",
            gtf=prepared / "annotation.parquet",
            max_mismatches=1,
            structure=StructureConfig(enabled=False),
            output_dir=tmp_path / "results",
        )
    )
    assert result.candidate_count == 6


def test_generic_transcriptome_preparation_keeps_provider_metadata(tmp_path: Path) -> None:
    prepared = tmp_path / "prepared-generic"

    manifest_path = prepare_transcriptome(
        DATA / "synthetic.fa",
        DATA / "synthetic.gtf",
        prepared,
        species="Mus musculus",
        assembly="GRCm39",
        provider="GENCODE",
        release="M39",
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["species"] == "Mus musculus"
    assert manifest["assembly"] == "GRCm39"
    assert manifest["release"] == "M39"


def test_zero_hit_run_still_writes_typed_empty_outputs(tmp_path: Path) -> None:
    result = run_scan(
        ScanConfig(
            query="TTTTTTTTTTTT",
            fasta=DATA / "synthetic.fa",
            gtf=DATA / "synthetic.gtf",
            max_mismatches=0,
            structure=StructureConfig(enabled=False),
            output_dir=tmp_path,
        )
    )
    hits = pd.read_parquet(result.output_dir / "all_transcript_hits.parquet")
    assert hits.empty
    assert {"rank", "risk_score", "transcript_id", "genomic_blocks"} <= set(hits.columns)
    assert pd.read_csv(result.output_dir / "brain_candidate_panel.csv").empty
