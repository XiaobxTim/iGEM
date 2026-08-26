from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import pandas as pd
from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from plotly.offline import get_plotlyjs

from pufscan.config import ScanConfig, StructureConfig
from pufscan.consequence import parse_editing_window
from pufscan.plotting import (
    accessibility_plot,
    mismatch_distribution,
    risk_distribution,
    transcript_schematic,
)
from pufscan.transcriptomes import TranscriptomeRegistry
from pufscan.web.jobs import JobStore
from pufscan.web.queries import get_candidate, query_candidates
from pufscan.web.uploads import (
    ALLOWED_FASTA_SUFFIXES,
    ALLOWED_GTF_SUFFIXES,
    _matching_suffix,
    discard_staging,
    new_staging_directory,
    prepare_custom_upload,
    save_upload,
)


def create_app(
    project_root: Path | None = None,
    database_path: Path | None = None,
    results_dir: Path | None = None,
    custom_data_dir: Path | None = None,
    max_upload_bytes: int = 512 * 1024 * 1024,
) -> FastAPI:
    root = project_root or Path.cwd()
    database = database_path or root / ".pufscan_web/atlas.sqlite3"
    results = results_dir or root / "results/web"
    custom_data = custom_data_dir or root / "data/custom"
    templates = Jinja2Templates(directory=Path(__file__).parent / "templates")

    app = FastAPI(title="PUF-OffTarget Atlas")
    app.mount("/static", StaticFiles(directory=Path(__file__).parent / "static"), name="static")
    app.state.project_root = root
    app.state.results_dir = results
    app.state.registry = TranscriptomeRegistry(database, root)
    app.state.jobs = JobStore(database)
    app.state.custom_data_dir = custom_data
    app.state.max_upload_bytes = max_upload_bytes

    @app.get("/", response_class=HTMLResponse)
    def home(request: Request) -> HTMLResponse:
        transcriptomes = [
            {"spec": spec, "availability": app.state.registry.availability(spec.id)}
            for spec in app.state.registry.list()
        ]
        return templates.TemplateResponse(
            request=request,
            name="home.html",
            context={"transcriptomes": transcriptomes, "error": None},
        )

    @app.get("/assets/plotly.js")
    def plotly_asset() -> Response:
        return Response(
            get_plotlyjs(),
            media_type="application/javascript",
            headers={"Cache-Control": "public, max-age=604800, immutable"},
        )

    @app.post("/runs", response_class=HTMLResponse)
    async def create_run(request: Request) -> Response:
        form = await request.form()
        transcriptome_id = str(form.get("transcriptome_id", ""))
        try:
            spec = app.state.registry.get(transcriptome_id)
            availability = app.state.registry.availability(transcriptome_id)
            if not availability.ready:
                raise ValueError("The selected transcriptome is not installed on this server")
            mode = str(form.get("mode", "binding_only"))
            window_raw = str(form.get("editing_window", ""))
            mismatch_raw = form.get("max_mismatches", "2")
            max_mismatches = int(mismatch_raw) if isinstance(mismatch_raw, str) else 2
            release_raw = spec.release.removeprefix("M")
            release_number = int(release_raw) if release_raw.isdigit() else 1
            config = ScanConfig(
                query=str(form.get("query", "")),
                fasta=spec.fasta_path,
                gtf=spec.annotation_path,
                expression=spec.expression_path if spec.expression_path and spec.expression_path.exists() else None,
                max_mismatches=max_mismatches,
                search_reverse_complement=form.get("search_reverse_complement") == "on",
                mode=mode,  # type: ignore[arg-type]
                editor=str(form.get("editor")) if mode == "editor_fusion" else None,  # type: ignore[arg-type]
                editing_window=parse_editing_window(window_raw) if mode == "editor_fusion" else None,
                structure=StructureConfig(enabled=form.get("structure_enabled") == "on"),
                output_dir=results,
                gencode_release=release_number,
                species=spec.species,
                genome_build=spec.assembly,
                annotation_provider=spec.provider,
                annotation_release=spec.release,
            )
        except Exception:
            transcriptomes = [
                {"spec": item, "availability": app.state.registry.availability(item.id)}
                for item in app.state.registry.list()
            ]
            return templates.TemplateResponse(
                request=request,
                name="home.html",
                context={
                    "transcriptomes": transcriptomes,
                    "error": "Enter 8–12 RNA bases (A, C, G, or U) and check the analysis settings.",
                },
                status_code=422,
            )
        record = app.state.jobs.create(config.model_dump(mode="json"), transcriptome_id)
        return RedirectResponse(url=f"/runs/{record.id}", status_code=303)

    @app.post("/transcriptomes/custom")
    async def upload_transcriptome(
        identifier: Annotated[str, Form()],
        display_name: Annotated[str, Form()],
        species: Annotated[str, Form()],
        assembly: Annotated[str, Form()],
        provider: Annotated[str, Form()],
        release: Annotated[str, Form()],
        fasta: Annotated[UploadFile, File()],
        gtf: Annotated[UploadFile, File()],
    ) -> Response:
        staging = new_staging_directory(custom_data)
        try:
            fasta_suffix = _matching_suffix(fasta.filename, ALLOWED_FASTA_SUFFIXES)
            gtf_suffix = _matching_suffix(gtf.filename, ALLOWED_GTF_SUFFIXES)
            source_dir = staging / "source"
            fasta_path = await save_upload(
                fasta, source_dir / f"transcripts{fasta_suffix}", max_upload_bytes
            )
            gtf_path = await save_upload(
                gtf, source_dir / f"annotation{gtf_suffix}", max_upload_bytes
            )
            spec = prepare_custom_upload(
                staging,
                custom_data / identifier,
                identifier,
                display_name,
                species,
                assembly,
                provider,
                release,
                fasta_path,
                gtf_path,
            )
            app.state.registry.add(spec)
        except OverflowError as error:
            discard_staging(staging)
            return HTMLResponse(str(error), status_code=413)
        except Exception as error:
            discard_staging(staging)
            return HTMLResponse(str(error), status_code=422)
        return RedirectResponse(url="/", status_code=303)

    @app.get("/runs/{run_id}", response_class=HTMLResponse)
    def run_page(request: Request, run_id: str) -> HTMLResponse:
        try:
            run = app.state.jobs.get(run_id)
        except KeyError:
            return HTMLResponse("Run not found", status_code=404)
        if run.status == "completed" and run.output_dir is not None:
            summary = json.loads((run.output_dir / "summary.json").read_text(encoding="utf-8"))
            hits = pd.read_parquet(run.output_dir / "all_transcript_hits.parquet")
            candidate_page = query_candidates(
                run.output_dir / "all_transcript_hits.parquet", page=1, page_size=25
            )
            return templates.TemplateResponse(
                request=request,
                name="results.html",
                context={
                    "run": run,
                    "summary": summary,
                    "candidate_page": candidate_page,
                    "mismatch_plot": mismatch_distribution(hits) if not hits.empty else None,
                    "risk_plot": risk_distribution(hits) if not hits.empty else None,
                },
            )
        return templates.TemplateResponse(request=request, name="run.html", context={"run": run})

    @app.get("/runs/{run_id}/progress", response_class=HTMLResponse)
    def run_progress(request: Request, run_id: str) -> HTMLResponse:
        try:
            run = app.state.jobs.get(run_id)
        except KeyError:
            return HTMLResponse("Run not found", status_code=404)
        return templates.TemplateResponse(
            request=request, name="partials/progress.html", context={"run": run}
        )

    @app.get("/runs/{run_id}/candidates", response_class=HTMLResponse)
    def candidates(
        request: Request,
        run_id: str,
        gene: str = "",
        transcript: str = "",
        region: str = "",
        mismatch: str = "",
        min_risk: float = 0,
        sort: str = "rank",
        direction: str = "asc",
        page: int = 1,
        page_size: int = 25,
    ) -> HTMLResponse:
        try:
            run = app.state.jobs.get(run_id)
            if run.status != "completed" or run.output_dir is None:
                raise KeyError(run_id)
            regions = tuple(value for value in region.split(",") if value)
            mismatches = tuple(int(value) for value in mismatch.split(",") if value)
            candidate_page = query_candidates(
                run.output_dir / "all_transcript_hits.parquet",
                gene=gene,
                transcript=transcript,
                regions=regions,
                mismatches=mismatches,
                min_risk=min_risk,
                sort=sort,
                direction=direction,
                page=page,
                page_size=page_size,
            )
        except (KeyError, ValueError):
            return HTMLResponse("Results not found", status_code=404)
        return templates.TemplateResponse(
            request=request,
            name="partials/candidate_table.html",
            context={"run": run, "candidate_page": candidate_page},
        )

    @app.get("/runs/{run_id}/candidates/{rank}", response_class=HTMLResponse)
    def candidate_detail(request: Request, run_id: str, rank: int) -> HTMLResponse:
        try:
            run = app.state.jobs.get(run_id)
            if run.status != "completed" or run.output_dir is None:
                raise KeyError(run_id)
            candidate = get_candidate(run.output_dir / "all_transcript_hits.parquet", rank)
        except KeyError:
            return HTMLResponse("Candidate not found", status_code=404)
        return templates.TemplateResponse(
            request=request,
            name="candidate_detail.html",
            context={
                "run": run,
                "candidate": candidate,
                "alignment": list(zip(candidate["query_rna"], candidate["matched_sequence_rna"], strict=True)),
                "schematic": transcript_schematic(pd.Series(candidate)),
                "accessibility": accessibility_plot(pd.Series(candidate)),
            },
        )

    download_allowlist = {
        "report.html",
        "top_hits.tsv",
        "all_transcript_hits.tsv.gz",
        "all_transcript_hits.parquet",
        "candidates.bed",
        "candidates.bed12",
        "potential_editing_events.tsv",
        "summary.json",
        "run_metadata.json",
    }

    @app.get("/runs/{run_id}/downloads/{artifact}")
    def download_result(run_id: str, artifact: str) -> Response:
        try:
            run = app.state.jobs.get(run_id)
            if artifact not in download_allowlist or run.output_dir is None:
                raise KeyError(artifact)
            path = run.output_dir / artifact
            if not path.is_file():
                raise KeyError(artifact)
        except KeyError:
            return HTMLResponse("Download not found", status_code=404)
        return FileResponse(path, filename=artifact)

    return app


app = create_app()
