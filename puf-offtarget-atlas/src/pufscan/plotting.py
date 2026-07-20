from __future__ import annotations

from typing import cast

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def mismatch_distribution(frame: pd.DataFrame) -> str:
    counts = frame["mismatch_count"].value_counts().reindex(range(4), fill_value=0).sort_index()
    figure = px.bar(x=counts.index, y=counts.values, labels={"x": "Mismatches", "y": "Candidates"})
    figure.update_layout(template="plotly_white", height=330, margin=dict(l=40, r=20, t=30, b=40))
    return cast(str, figure.to_html(full_html=False, include_plotlyjs=False))


def risk_distribution(frame: pd.DataFrame) -> str:
    figure = px.histogram(frame, x="risk_score", nbins=20, labels={"risk_score": "Risk priority score"})
    figure.update_layout(template="plotly_white", height=330, margin=dict(l=40, r=20, t=30, b=40))
    return cast(str, figure.to_html(full_html=False, include_plotlyjs=False))


def transcript_schematic(row: pd.Series) -> str:
    import json

    end = int(row["transcript_end"])
    start = int(row["transcript_start"])
    features = json.loads(row.get("transcript_features", "{}"))
    junctions = json.loads(row.get("transcript_junctions", "[]"))
    transcript_length = max([end, *[interval[1] for intervals in features.values() for interval in intervals]])
    figure = go.Figure()
    figure.add_trace(go.Scatter(x=[1, max(transcript_length, 2)], y=[0, 0], mode="lines", name="Transcript", line=dict(color="#65757c", width=7)))
    colors = {"UTR": "#9bbfc3", "five_prime_utr": "#9bbfc3", "three_prime_utr": "#9bbfc3", "CDS": "#087f8c"}
    for feature, color in colors.items():
        for feature_start, feature_end in features.get(feature, []):
            figure.add_trace(go.Scatter(x=[feature_start + 1, feature_end], y=[0, 0], mode="lines", name=feature, line=dict(color=color, width=18), showlegend=False))
    for junction in junctions:
        figure.add_vline(x=junction, line_dash="dot", line_color="#6d7478")
    figure.add_vrect(x0=start, x1=end, fillcolor="#d94841", opacity=0.45, line_width=0, annotation_text="PUF motif")
    positions = json.loads(row.get("potential_edit_positions", "[]"))
    if positions:
        figure.add_trace(go.Scatter(x=positions, y=[0.08] * len(positions), mode="markers", name="Potential editable base", marker=dict(color="#f0a202", size=9)))
    figure.update_yaxes(visible=False)
    figure.update_layout(template="plotly_white", height=220, margin=dict(l=30, r=20, t=30, b=30))
    return cast(str, figure.to_html(full_html=False, include_plotlyjs=False))


def tissue_heatmap(frame: pd.DataFrame, top_n: int = 20) -> str:
    import json

    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in frame.sort_values("risk_score", ascending=False).to_dict(orient="records"):
        gene = str(row.get("gene_name") or row.get("gene_id"))
        if gene in seen:
            continue
        values = json.loads(row.get("all_tissue_tpm") or "null")
        if values:
            rows.append({"gene": gene, **values})
            seen.add(gene)
        if len(rows) >= top_n:
            break
    if not rows:
        return "<p>Not available</p>"
    heatmap = pd.DataFrame(rows).set_index("gene")
    figure = px.imshow(heatmap, aspect="auto", color_continuous_scale="Viridis", labels={"color": "TPM"})
    figure.update_layout(template="plotly_white", height=max(300, 28 * len(heatmap)), margin=dict(l=70, r=20, t=30, b=80))
    return cast(str, figure.to_html(full_html=False, include_plotlyjs=False))


def accessibility_plot(row: pd.Series) -> str:
    import json

    profile = json.loads(row.get("accessibility_profile") or "[]")
    if not profile or not any(value is not None for value in profile):
        return "<p>Not available</p>"
    window_start = int(row.get("structure_window_transcript_start") or 1)
    positions = list(range(window_start, window_start + len(profile)))
    figure = px.line(x=positions, y=profile, labels={"x": "Transcript position", "y": "Unpaired probability"})
    figure.add_vrect(x0=row["transcript_start"], x1=row["transcript_end"], fillcolor="#d94841", opacity=0.2, line_width=0)
    figure.update_yaxes(range=[0, 1])
    figure.update_layout(template="plotly_white", height=330, margin=dict(l=50, r=20, t=30, b=40))
    return cast(str, figure.to_html(full_html=False, include_plotlyjs=False))
