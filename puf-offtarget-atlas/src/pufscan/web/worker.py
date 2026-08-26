from __future__ import annotations

import time

from pufscan.config import ScanConfig
from pufscan.pipeline import run_scan
from pufscan.web.jobs import JobStore, RunStatus


def _status_for_progress(message: str, fraction: float) -> RunStatus:
    lowered = message.lower()
    if "search" in lowered:
        return "searching"
    if "annotation" in lowered or "annotating" in lowered:
        return "annotating"
    if "writing" in lowered or fraction >= 0.75:
        return "reporting"
    return "scoring"


def run_one_job(store: JobStore) -> bool:
    record = store.claim_next()
    if record is None:
        return False

    def update(message: str, fraction: float) -> None:
        store.update_progress(
            record.id, message, fraction, _status_for_progress(message, fraction)
        )

    try:
        result = run_scan(ScanConfig.model_validate(record.config), progress=update)
        store.complete(record.id, result.output_dir)
    except Exception as error:
        store.fail(record.id, str(error))
    return True


def run_worker(store: JobStore, poll_seconds: float = 1.0) -> None:
    while True:
        if not run_one_job(store):
            time.sleep(poll_seconds)


__all__ = ["run_one_job", "run_worker"]
