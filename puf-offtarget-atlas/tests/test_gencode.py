from pathlib import Path

from pufscan.gencode import parse_gtf, strip_version

DATA = Path(__file__).parent / "data"


def test_strip_version_preserves_stable_identifier() -> None:
    assert strip_version("ENST000001.12") == "ENST000001"


def test_parse_gtf_retains_required_fields_and_repeated_tags() -> None:
    records = list(parse_gtf(DATA / "synthetic.gtf"))
    transcript = next(row for row in records if row.feature == "transcript" and row.transcript_id == "ENST000001.1")
    assert transcript.gene_id_without_version == "ENSG000001"
    assert transcript.transcript_id_without_version == "ENST000001"
    assert transcript.tags == ("basic", "MANE_Select")
    cds = next(row for row in records if row.feature == "CDS" and row.transcript_id == "ENST000001.1")
    assert cds.phase == 0

