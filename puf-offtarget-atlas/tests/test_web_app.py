from pathlib import Path

from fastapi.testclient import TestClient

from pufscan.transcriptomes import TranscriptomeSpec
from pufscan.web.app import create_app
from pufscan.web.worker import run_one_job

DATA = Path(__file__).parent / "data"


def make_client(tmp_path: Path) -> TestClient:
    app = create_app(
        project_root=Path.cwd(),
        database_path=tmp_path / "atlas.sqlite3",
        results_dir=tmp_path / "results",
    )
    app.state.registry.add(
        TranscriptomeSpec(
            id="synthetic-demo",
            display_name="Synthetic demo",
            species="Testus organismus",
            assembly="Demo1",
            provider="Test",
            release="1",
            fasta_path=DATA / "synthetic.fa",
            annotation_path=DATA / "synthetic.gtf",
            source="custom",
        )
    )
    return TestClient(app)


def test_home_lists_available_transcriptomes(tmp_path: Path) -> None:
    client = make_client(tmp_path)

    response = client.get("/")

    assert response.status_code == 200
    assert "PUF-OffTarget Atlas" in response.text
    assert "Synthetic demo" in response.text
    assert "Advanced analysis settings" in response.text
    assert "Project Wiki" in response.text
    assert "unpkg.com" not in response.text
    assert 'class="base base-a"' in response.text
    assert client.get("/static/atlas.css").status_code == 200


def test_submit_scan_creates_queued_job(tmp_path: Path) -> None:
    client = make_client(tmp_path)

    response = client.post(
        "/runs",
        data={
            "query": "AACGUCUAUA",
            "transcriptome_id": "synthetic-demo",
            "max_mismatches": "1",
            "mode": "binding_only",
        },
        follow_redirects=False,
    )

    assert response.status_code == 303
    run_page = client.get(response.headers["location"])
    assert run_page.status_code == 200
    assert "Queued" in run_page.text


def test_submit_scan_returns_form_error_for_invalid_query(tmp_path: Path) -> None:
    client = make_client(tmp_path)

    response = client.post(
        "/runs",
        data={
            "query": "BAD",
            "transcriptome_id": "synthetic-demo",
            "max_mismatches": "1",
            "mode": "binding_only",
        },
    )

    assert response.status_code == 422
    assert "8–12 RNA bases" in response.text


def test_completed_run_exposes_results_and_candidate_partial(tmp_path: Path) -> None:
    client = make_client(tmp_path)
    response = client.post(
        "/runs",
        data={
            "query": "AACGUCUAUA",
            "transcriptome_id": "synthetic-demo",
            "max_mismatches": "1",
            "mode": "binding_only",
        },
        follow_redirects=False,
    )
    run_id = response.headers["location"].rsplit("/", 1)[-1]
    assert run_one_job(client.app.state.jobs) is True

    result = client.get(f"/runs/{run_id}")
    candidates = client.get(f"/runs/{run_id}/candidates?gene=GENE1")
    detail = client.get(f"/runs/{run_id}/candidates/1")

    assert "Candidate sites" in result.text
    assert "Use in Brain Delivery Model" in result.text
    assert f'/runs/{run_id}/downloads/brain_candidate_panel.csv' in result.text
    assert "/assets/plotly.js" in result.text
    assert "Filter candidate sites" in result.text
    assert candidates.status_code == 200
    assert "GENE1" in candidates.text
    assert detail.status_code == 200
    assert "Sequence alignment" in detail.text

    report_download = client.get(f"/runs/{run_id}/downloads/report.html")
    panel_download = client.get(f"/runs/{run_id}/downloads/brain_candidate_panel.csv")
    panel_metadata = client.get(
        f"/runs/{run_id}/downloads/brain_candidate_panel.metadata.json"
    )
    missing_download = client.get(f"/runs/{run_id}/downloads/secrets.txt")
    assert report_download.status_code == 200
    assert panel_download.status_code == 200
    assert panel_metadata.status_code == 200
    assert missing_download.status_code == 404


def test_progress_partial_polls_while_run_is_queued(tmp_path: Path) -> None:
    client = make_client(tmp_path)
    response = client.post(
        "/runs",
        data={
            "query": "AACGUCUAUA",
            "transcriptome_id": "synthetic-demo",
            "max_mismatches": "1",
            "mode": "binding_only",
        },
        follow_redirects=False,
    )
    run_id = response.headers["location"].rsplit("/", 1)[-1]

    progress = client.get(f"/runs/{run_id}/progress")

    assert progress.status_code == 200
    assert 'hx-trigger="every 2s"' in progress.text
    assert "Queued" in progress.text
