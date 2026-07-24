from __future__ import annotations

from models.editing.offtarget_panel import summarize_offtarget_panel


def test_offtarget_panel_summary_uses_accessibility_weighting():
    rows = [
        {
            "site_id": "a",
            "gene": "A",
            "initial_pool": 1.0,
            "binding_score": 0.5,
            "accessibility": 1.0,
            "context_score": 0.5,
        },
        {
            "site_id": "b",
            "gene": "B",
            "initial_pool": 1.0,
            "binding_score": 1.0,
            "accessibility": 0.0,
            "context_score": 1.0,
        },
    ]

    summary = summarize_offtarget_panel(rows)

    assert summary["effective_pool"] == 1.0
    assert summary["relative_binding"] == 0.5
    assert summary["relative_catalysis"] == 0.25
