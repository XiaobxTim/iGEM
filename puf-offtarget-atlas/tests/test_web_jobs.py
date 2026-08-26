from pathlib import Path

from pufscan.config import ScanConfig, StructureConfig
from pufscan.web.jobs import JobStore
from pufscan.web.worker import run_one_job

DATA = Path(__file__).parent / "data"


def test_job_store_claims_oldest_queued_run_and_persists_progress(tmp_path: Path) -> None:
    store = JobStore(tmp_path / "atlas.sqlite3")
    first = store.create({"query": "AACGTCTATA"}, "human-gencode50")
    store.create({"query": "TTTTTTTT"}, "mouse-gencode-m39")

    claimed = store.claim_next()
    assert claimed is not None
    assert claimed.id == first.id
    assert claimed.status == "preparing"

    store.update_progress(first.id, "Searching transcriptome", 0.2, "searching")
    updated = store.get(first.id)
    assert updated.progress == 0.2
    assert updated.stage == "Searching transcriptome"
    assert updated.status == "searching"


def test_job_store_records_completion_and_failure(tmp_path: Path) -> None:
    store = JobStore(tmp_path / "atlas.sqlite3")
    completed = store.create({"query": "AACGTCTATA"}, "human-gencode50")
    failed = store.create({"query": "TTTTTTTT"}, "human-gencode50")

    store.complete(completed.id, tmp_path / "result")
    store.fail(failed.id, "RNAplfold exploded")

    assert store.get(completed.id).status == "completed"
    assert store.get(completed.id).output_dir == tmp_path / "result"
    assert store.get(failed.id).status == "failed"
    assert store.get(failed.id).error == "RNAplfold exploded"


def test_worker_runs_scan_and_records_output(tmp_path: Path) -> None:
    store = JobStore(tmp_path / "atlas.sqlite3")
    config = ScanConfig(
        query="AACGUCUAUA",
        fasta=DATA / "synthetic.fa",
        gtf=DATA / "synthetic.gtf",
        max_mismatches=0,
        structure=StructureConfig(enabled=False),
        output_dir=tmp_path / "results",
    )
    queued = store.create(config.model_dump(mode="json"), "synthetic-demo")

    assert run_one_job(store) is True

    finished = store.get(queued.id)
    assert finished.status == "completed"
    assert finished.output_dir is not None
    assert (finished.output_dir / "summary.json").exists()
