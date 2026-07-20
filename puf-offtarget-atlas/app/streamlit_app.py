from __future__ import annotations

import io
import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from pufscan.config import ScanConfig, StructureConfig
from pufscan.consequence import parse_editing_window
from pufscan.pipeline import run_scan

st.set_page_config(page_title="PUF-OffTarget Atlas", layout="wide")
st.markdown(
    """
    <style>
    html, body, [class*="css"] { letter-spacing: 0 !important; }
    .block-container { max-width: 1500px; padding-top: 1.4rem; }
    h1 { font-family: Georgia, serif; font-size: 2.25rem !important; }
    [data-testid="stMetric"] { border-top: 3px solid #087f8c; padding-top: .5rem; }
    div[data-testid="stDataFrame"] { border: 1px solid #d9e0e3; }
    </style>
    """,
    unsafe_allow_html=True,
)
st.title("PUF-OffTarget Atlas")
st.caption("Potential transcriptome binding candidates · heuristic prioritization score, not a calibrated probability")

with st.sidebar:
    st.header("Analysis")
    query = st.text_input("PUF target RNA sequence", "AACGUCUAUA")
    mismatches = st.slider("Maximum mismatches", 0, 3, 2)
    fasta = st.text_input("GENCODE transcript FASTA", "tests/data/synthetic.fa")
    gtf = st.text_input("GENCODE annotation GTF", "tests/data/synthetic.gtf")
    expression = st.text_input("GTEx expression file", "")
    tissues = st.text_input("Target tissue", "")
    mode = st.selectbox("Analysis mode", ["binding_only", "editor_fusion"])
    editor = None
    editing_window = None
    if mode == "editor_fusion":
        editor = st.selectbox("Editor type", ["APOBEC_C2U", "ADAR_A2I"])
        editing_window = st.text_input("Editing window", "-15:10")
    structure_enabled = st.toggle("Structure analysis", value=False)
    structure_top_n = st.number_input("Top N structure candidates", min_value=0, value=5000)
    position_weights = st.text_input("Position-weight file", "")
    output = st.text_input("Output directory", "results")
    run_clicked = st.button("Run analysis", type="primary", use_container_width=True)

if run_clicked:
    status = st.status("Preparing analysis", expanded=True)
    progress_bar = st.progress(0.0)

    def update_progress(message: str, fraction: float) -> None:
        status.write(message)
        progress_bar.progress(fraction)

    try:
        config = ScanConfig(
            query=query,
            fasta=Path(fasta),
            gtf=Path(gtf),
            expression=Path(expression) if expression else None,
            target_tissues=tuple(value.strip() for value in tissues.split(",") if value.strip()),
            max_mismatches=mismatches,
            mode=mode,
            editor=editor,
            editing_window=parse_editing_window(editing_window) if editing_window else None,
            structure=StructureConfig(enabled=structure_enabled, top_n=int(structure_top_n)),
            position_weights=Path(position_weights) if position_weights else None,
            output_dir=Path(output),
        )
        result = run_scan(config, update_progress)
        st.session_state["run_dir"] = str(result.output_dir)
        status.update(label="Analysis complete", state="complete", expanded=False)
    except Exception as error:
        status.update(label="Analysis failed", state="error")
        st.error(str(error))

run_dir_value = st.session_state.get("run_dir")
if run_dir_value:
    run_dir = Path(run_dir_value)
    hits = pd.read_parquet(run_dir / "all_transcript_hits.parquet")
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    metric_columns = st.columns(4)
    for column, key in zip(metric_columns, ["candidate_sites", "candidate_genes", "candidate_transcripts", "unique_genomic_loci"], strict=True):
        column.metric(key.replace("_", " ").title(), summary[key])
    st.subheader("Candidates")
    filters = st.columns(5)
    genes = filters[0].multiselect("Gene", sorted(hits["gene_name"].dropna().unique()))
    transcripts = filters[1].multiselect("Transcript", sorted(hits["transcript_id"].unique()))
    regions = filters[2].multiselect("Region", sorted(hits["transcript_region"].unique()))
    mismatch_filter = filters[3].multiselect("Mismatch", sorted(hits["mismatch_count"].unique()))
    minimum_risk = filters[4].number_input("Minimum risk score", 0.0, 100.0, 0.0)
    filtered = hits.copy()
    if genes:
        filtered = filtered[filtered["gene_name"].isin(genes)]
    if transcripts:
        filtered = filtered[filtered["transcript_id"].isin(transcripts)]
    if regions:
        filtered = filtered[filtered["transcript_region"].isin(regions)]
    if mismatch_filter:
        filtered = filtered[filtered["mismatch_count"].isin(mismatch_filter)]
    filtered = filtered[filtered["risk_score"] >= minimum_risk]
    event = st.dataframe(filtered, hide_index=True, use_container_width=True, on_select="rerun", selection_mode="single-row")
    download_columns = st.columns(4)
    download_columns[0].download_button("CSV", filtered.to_csv(index=False), "puf_candidates.csv", "text/csv")
    download_columns[1].download_button("TSV", filtered.to_csv(index=False, sep="\t"), "puf_candidates.tsv", "text/tab-separated-values")
    parquet_buffer = io.BytesIO()
    filtered.to_parquet(parquet_buffer, index=False)
    download_columns[2].download_button("Parquet", parquet_buffer.getvalue(), "puf_candidates.parquet", "application/octet-stream")
    download_columns[3].download_button("HTML report", (run_dir / "report.html").read_bytes(), "puf_report.html", "text/html")
    selected_indices = event.selection.rows
    if selected_indices:
        row = filtered.iloc[selected_indices[0]]
        st.subheader(f"{row['gene_name']} · {row['transcript_id']}")
        left, right = st.columns([1, 2])
        left.code(row["match_pattern"], language=None)
        left.write(row["consequence_summary"])
        expression_values = json.loads(row.get("all_tissue_tpm") or "null")
        if expression_values:
            right.plotly_chart(px.bar(x=list(expression_values), y=list(expression_values.values()), labels={"x": "Tissue", "y": "TPM"}), use_container_width=True)
        profile = json.loads(row.get("accessibility_profile") or "[]")
        if profile and any(value is not None for value in profile):
            st.plotly_chart(px.line(y=profile, labels={"index": "Local position", "value": "Unpaired probability"}), use_container_width=True)
        events_path = run_dir / "potential_editing_events.tsv"
        editing_events = pd.read_csv(events_path, sep="\t")
        if not editing_events.empty:
            st.dataframe(editing_events[editing_events["transcript_id"] == row["transcript_id"]], hide_index=True, use_container_width=True)
