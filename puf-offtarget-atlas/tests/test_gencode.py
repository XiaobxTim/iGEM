from pathlib import Path

from pufscan.gencode import gencode_filenames, gencode_release_url, parse_gtf, strip_version

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


def test_mouse_release_uses_mouse_ftp_layout() -> None:
    assert gencode_release_url("M39", "mouse").endswith("Gencode_mouse/release_M39")
    assert gencode_filenames("M39", "mouse")["transcript_fasta"] == "gencode.vM39.transcripts.fa.gz"
