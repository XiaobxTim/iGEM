from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from plotly.offline import get_plotlyjs
from starlette.datastructures import UploadFile

from brain_app.service import get_designs, parse_candidate_panel, run_model


def create_app() -> FastAPI:
    package_dir = Path(__file__).resolve().parent
    templates = Jinja2Templates(directory=package_dir / "templates")
    app = FastAPI(title="Brain Delivery Model", version="4.0")
    app.mount("/static", StaticFiles(directory=package_dir / "static"), name="static")

    @app.get("/", response_class=HTMLResponse)
    def home(request: Request) -> HTMLResponse:
        return templates.TemplateResponse(
            request=request,
            name="index.html",
            context={
                "designs": get_designs(),
                "wiki_url": os.getenv("IGEM_WIKI_URL", "http://127.0.0.1:5173"),
            },
        )

    @app.get("/assets/plotly.js")
    def plotly_asset() -> Response:
        return Response(
            get_plotlyjs(),
            media_type="application/javascript",
            headers={"Cache-Control": "public, max-age=604800, immutable"},
        )

    @app.get("/healthz")
    def health() -> dict[str, str]:
        return {"status": "ok", "model": "brain-delivery-v4"}

    @app.post("/api/simulate")
    async def simulate(request: Request) -> Response:
        try:
            form = await request.form()
            upload = form.get("candidate_panel")
            panel_rows = None
            if isinstance(upload, UploadFile) and upload.filename:
                panel_rows = parse_candidate_panel(await upload.read())
            payload = run_model(
                mode=str(form.get("mode", "single")),
                design_id=str(form.get("design_id", "")),
                route=str(form.get("route", "footpad")),
                dose=float(str(form.get("dose", "1.0"))),
                duration=float(str(form.get("duration", "48"))),
                panel_rows=panel_rows,
            )
            return JSONResponse(payload)
        except (TypeError, ValueError) as error:
            return JSONResponse({"detail": str(error)}, status_code=422)
        except RuntimeError as error:
            return JSONResponse({"detail": str(error)}, status_code=500)

    return app


app = create_app()
