from __future__ import annotations

import io

from fastapi.testclient import TestClient

from brain_app.app import create_app
from brain_app.service import MAX_PANEL_BYTES, parse_candidate_panel, run_model


def test_home_exposes_model_inputs_and_disclaimer() -> None:
    client = TestClient(create_app())

    response = client.get("/")

    assert response.status_code == 200
    assert "Brain Delivery Model" in response.text
    assert "Single simulation" in response.text
    assert "Dose optimization" in response.text
    assert "not clinically calibrated" in response.text
    assert client.get("/assets/plotly.js").status_code == 200


def test_candidate_panel_parser_enforces_contract() -> None:
    valid = (
        "site_id,gene,initial_pool,binding_score,accessibility,context_score,"
        "validation_priority,notes\n"
        "puf-1,APOE,0.8,0.9,0.7,0.6,High priority,synthetic\n"
    ).encode()

    rows = parse_candidate_panel(valid)

    assert rows[0]["site_id"] == "puf-1"
    assert rows[0]["binding_score"] == 0.9


def test_candidate_panel_parser_rejects_invalid_scores_and_large_uploads() -> None:
    invalid = (
        "site_id,gene,initial_pool,binding_score,accessibility,context_score,"
        "validation_priority\n"
        "bad,APOE,0.8,1.2,0.7,0.6,High\n"
    ).encode()

    try:
        parse_candidate_panel(invalid)
    except ValueError as error:
        assert "binding_score" in str(error)
    else:
        raise AssertionError("invalid score was accepted")

    try:
        parse_candidate_panel(b"x" * (MAX_PANEL_BYTES + 1))
    except ValueError as error:
        assert "5 MB" in str(error)
    else:
        raise AssertionError("oversized panel was accepted")


def test_single_simulation_api_applies_design_and_candidate_panel() -> None:
    client = TestClient(create_app())
    panel = (
        "site_id,gene,initial_pool,binding_score,accessibility,context_score,"
        "validation_priority,notes\n"
        "puf-1,GENE1,0.8,0.9,0.7,0.6,High priority,synthetic\n"
    ).encode()

    response = client.post(
        "/api/simulate",
        data={
            "mode": "single",
            "design_id": "10R-design-4",
            "route": "footpad",
            "dose": "1.0",
            "duration": "24",
        },
        files={"candidate_panel": ("brain_candidate_panel.csv", io.BytesIO(panel), "text/csv")},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["mode"] == "single"
    assert payload["inputs"]["design_id"] == "10R-design-4"
    assert payload["panel_summary"]["n_sites"] == 1.0
    assert 0.0 <= payload["metrics"]["apoe3_like_fraction_final"] <= 1.0
    assert 0.0 <= payload["metrics"]["apoe2_like_fraction_final"] <= 1.0
    assert len(payload["series"]["time"]) == len(payload["series"]["pbrain"])


def test_dose_optimization_api_returns_twelve_doses() -> None:
    client = TestClient(create_app())

    response = client.post(
        "/api/simulate",
        data={
            "mode": "optimize",
            "design_id": "10R-design-4",
            "route": "im",
            "dose": "1.0",
            "duration": "24",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["mode"] == "optimize"
    assert len(payload["dose_scan"]) == 12
    assert payload["dose_scan"][0]["dose"] == 0.1


def test_invalid_request_returns_explainable_error() -> None:
    client = TestClient(create_app())

    response = client.post(
        "/api/simulate",
        data={
            "mode": "single",
            "design_id": "not-a-design",
            "route": "iv",
            "dose": "99",
            "duration": "24",
        },
    )

    assert response.status_code == 422
    assert "dose" in response.json()["detail"].lower()


def test_empty_uploaded_panel_retains_conservative_base_prior() -> None:
    inputs = {
        "mode": "single",
        "design_id": "10R-design-4",
        "route": "footpad",
        "dose": 1.0,
        "duration": 24.0,
    }

    baseline = run_model(**inputs, panel_rows=None)
    empty_panel = run_model(**inputs, panel_rows=[])

    assert empty_panel["metrics"] == baseline["metrics"]
    assert empty_panel["panel_summary"]["provided"] is True
    assert "retained" in empty_panel["panel_summary"]["warning"].lower()


def test_web_metrics_match_direct_service_call() -> None:
    client = TestClient(create_app())
    form = {
        "mode": "single",
        "design_id": "10R-design-4",
        "route": "iv",
        "dose": "0.8",
        "duration": "24",
    }

    response = client.post("/api/simulate", data=form)
    direct = run_model(
        mode="single",
        design_id="10R-design-4",
        route="iv",
        dose=0.8,
        duration=24.0,
    )

    assert response.status_code == 200
    assert response.json()["metrics"] == direct["metrics"]


def test_web_rejects_oversized_panel_before_parsing() -> None:
    client = TestClient(create_app())

    response = client.post(
        "/api/simulate",
        data={
            "mode": "single",
            "design_id": "10R-design-4",
            "route": "footpad",
            "dose": "1.0",
            "duration": "24",
        },
        files={"candidate_panel": ("too-large.csv", b"x" * (MAX_PANEL_BYTES + 1), "text/csv")},
    )

    assert response.status_code == 422
    assert "5 MB" in response.json()["detail"]
