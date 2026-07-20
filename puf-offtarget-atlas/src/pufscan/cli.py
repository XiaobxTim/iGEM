from __future__ import annotations

import csv
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from pufscan.config import ScanConfig, StructureConfig, load_yaml
from pufscan.consequence import parse_editing_window
from pufscan.expression import prepare_expression
from pufscan.gencode import download_gencode
from pufscan.index import prepare_gencode
from pufscan.pipeline import run_scan
from pufscan.report import generate_report

app = typer.Typer(help="Prioritize potential PUF RNA-binding candidates in the human transcriptome.", no_args_is_help=True)
LOGGER = logging.getLogger(__name__)


@app.callback()
def main(verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False) -> None:
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO, format="%(levelname)s %(message)s")


@app.command("download-gencode")
def download_gencode_command(
    release: Annotated[int, typer.Option()] = 50,
    output: Annotated[Path, typer.Option()] = Path("data/gencode_v50"),
    all_regions: Annotated[bool, typer.Option("--all-regions/--reference-only")] = True,
) -> None:
    manifest = download_gencode(release, output, all_regions)
    typer.echo(f"Download manifest: {manifest}")


@app.command("prepare-gencode")
def prepare_gencode_command(
    fasta: Annotated[Path, typer.Option("--fasta")],
    gtf: Annotated[Path, typer.Option("--gtf")],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    manifest = prepare_gencode(fasta, gtf, output)
    typer.echo(f"Prepared manifest: {manifest}")


@app.command("prepare-gtex")
def prepare_gtex_command(
    input_path: Annotated[Path, typer.Option("--input")],
    output: Annotated[Path, typer.Option("--output")],
    level: Annotated[str, typer.Option()] = "gene",
) -> None:
    id_column = "gene_id" if level == "gene" else "transcript_id"
    prepare_expression(input_path, output, id_column=id_column)
    typer.echo(f"Prepared expression matrix: {output}")


@app.command("scan")
def scan_command(
    query: Annotated[str | None, typer.Option("--query")] = None,
    gencode_fasta: Annotated[Path | None, typer.Option("--gencode-fasta")] = None,
    gencode_gtf: Annotated[Path | None, typer.Option("--gencode-gtf")] = None,
    config_file: Annotated[Path | None, typer.Option("--config")] = None,
    max_mismatches: Annotated[int | None, typer.Option("--max-mismatches")] = None,
    gtex_expression: Annotated[Path | None, typer.Option("--gtex-expression")] = None,
    target_tissue: Annotated[list[str] | None, typer.Option("--target-tissue")] = None,
    mode: Annotated[str | None, typer.Option()] = None,
    editor: Annotated[str | None, typer.Option()] = None,
    editing_window: Annotated[str | None, typer.Option("--editing-window")] = None,
    structure: Annotated[bool, typer.Option("--structure/--no-structure")] = True,
    structure_top_n: Annotated[int, typer.Option("--structure-top-n")] = 5000,
    search_reverse_complement: Annotated[bool, typer.Option("--search-reverse-complement")] = False,
    position_weight_file: Annotated[Path | None, typer.Option("--position-weight-file")] = None,
    substitution_matrix: Annotated[Path | None, typer.Option("--substitution-matrix")] = None,
    gene_metadata: Annotated[Path | None, typer.Option("--gene-metadata")] = None,
    ensembl_rest: Annotated[bool, typer.Option("--ensembl-rest")] = False,
    output: Annotated[Path, typer.Option("--output")] = Path("results"),
    threads: Annotated[int, typer.Option()] = 1,
) -> None:
    file_config = load_yaml(config_file) if config_file else {}
    query = query or file_config.get("query")
    gencode_fasta = gencode_fasta or (
        Path(file_config["fasta"]) if file_config.get("fasta") else None
    )
    gencode_gtf = gencode_gtf or (Path(file_config["gtf"]) if file_config.get("gtf") else None)
    if query is None or gencode_fasta is None or gencode_gtf is None:
        raise typer.BadParameter("query, gencode FASTA, and gencode GTF are required via CLI or YAML")
    max_mismatches = (
        max_mismatches if max_mismatches is not None else int(file_config.get("max_mismatches", 2))
    )
    mode = mode or str(file_config.get("mode", "binding_only"))
    if gtex_expression is None and file_config.get("expression"):
        gtex_expression = Path(file_config["expression"])
    if not target_tissue and file_config.get("target_tissues"):
        target_tissue = list(file_config["target_tissues"])
    parsed_window = parse_editing_window(editing_window) if editing_window is not None else None
    config = ScanConfig(
        query=query,
        fasta=gencode_fasta,
        gtf=gencode_gtf,
        expression=gtex_expression,
        target_tissues=tuple(target_tissue or ()),
        max_mismatches=max_mismatches,
        search_reverse_complement=search_reverse_complement,
        mode=mode,  # type: ignore[arg-type]
        editor=editor,  # type: ignore[arg-type]
        editing_window=parsed_window,
        structure=StructureConfig(enabled=structure, top_n=structure_top_n),
        position_weights=position_weight_file,
        substitution_matrix=substitution_matrix,
        gene_metadata=gene_metadata,
        ensembl_rest=ensembl_rest,
        output_dir=output,
        threads=threads,
    )
    result = run_scan(config)
    typer.echo(f"Output directory: {result.output_dir}")
    typer.echo(f"candidate sites: {result.candidate_count}")


@app.command("report")
def report_command(run_directory: Annotated[Path, typer.Option("--run-directory")]) -> None:
    typer.echo(f"Report: {generate_report(run_directory)}")


@app.command("doctor")
def doctor_command(
    gencode_fasta: Annotated[Path | None, typer.Option("--gencode-fasta")] = None,
    gencode_gtf: Annotated[Path | None, typer.Option("--gencode-gtf")] = None,
    gtex_expression: Annotated[Path | None, typer.Option("--gtex-expression")] = None,
    output: Annotated[Path, typer.Option("--output")] = Path("results"),
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    checks = {
        "Python version": sys.version.split()[0],
        "GENCODE FASTA": "not specified" if gencode_fasta is None else ("OK" if gencode_fasta.exists() else "MISSING"),
        "FASTA index": "not specified" if gencode_fasta is None else ("OK" if Path(str(gencode_fasta) + ".fai").exists() else "not built"),
        "GENCODE GTF/index": "not specified" if gencode_gtf is None else ("OK" if gencode_gtf.exists() else "MISSING"),
        "GTEx file": "optional/not specified" if gtex_expression is None else ("OK" if gtex_expression.exists() else "MISSING"),
        "RNAplfold": shutil.which("RNAplfold") or "optional/not installed",
        "Output writable": "OK" if os.access(output, os.W_OK) else "NOT WRITABLE",
    }
    for name, value in checks.items():
        typer.echo(f"{name}: {value}")
    required_missing = any(value == "MISSING" for key, value in checks.items() if key.startswith("GENCODE"))
    if required_missing or checks["Output writable"] != "OK":
        raise typer.Exit(1)


@app.command("compare-designs")
def compare_designs_command(
    queries: Annotated[Path, typer.Option("--queries")],
    gencode_fasta: Annotated[Path, typer.Option("--gencode-fasta")],
    gencode_gtf: Annotated[Path, typer.Option("--gencode-gtf")],
    output: Annotated[Path, typer.Option("--output")],
) -> None:
    query_frame = pd.read_csv(queries)
    if "query" not in query_frame.columns:
        raise typer.BadParameter("Query CSV must contain a 'query' column")
    rows: list[dict[str, object]] = []
    for query in query_frame["query"].astype(str):
        result = run_scan(
            ScanConfig(query=query, fasta=gencode_fasta, gtf=gencode_gtf, structure=StructureConfig(enabled=False), output_dir=output / "runs")
        )
        hits = pd.read_parquet(result.output_dir / "all_transcript_hits.parquet")
        top = hits.head(10)["risk_score"] if not hits.empty else pd.Series(dtype=float)
        rows.append(
            {
                "query": query,
                "exact_transcript_hits": int((hits["mismatch_count"] == 0).sum()),
                "one_mismatch_hits": int((hits["mismatch_count"] == 1).sum()),
                "two_mismatch_hits": int((hits["mismatch_count"] == 2).sum()),
                "high_priority_genes": int(hits.loc[hits["risk_score"] >= 50, "gene_id"].nunique()),
                "high_priority_cds_hits": int(((hits["risk_score"] >= 50) & (hits["transcript_region"] == "CDS")).sum()),
                "target_tissue_exposure": hits["target_tissue_tpm"].max(),
                "maximum_off_target_risk": hits["risk_score"].max(),
                "mean_top_10_risk": top.mean() if len(top) else None,
            }
        )
    output.mkdir(parents=True, exist_ok=True)
    result_frame = pd.DataFrame(rows).sort_values("maximum_off_target_risk", na_position="last")
    result_frame["recommended_rank"] = range(1, len(result_frame) + 1)
    result_frame.to_csv(output / "design_comparison.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


if __name__ == "__main__":
    app()
