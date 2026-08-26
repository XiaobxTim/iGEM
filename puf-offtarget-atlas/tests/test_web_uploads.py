from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from pufscan.web.app import create_app
from pufscan.web.uploads import transcript_id_overlap

DATA = Path(__file__).parent / "data"


def test_transcript_overlap_rejects_incompatible_annotation(tmp_path: Path) -> None:
    fasta = tmp_path / "transcripts.fa"
    gtf = tmp_path / "annotation.gtf"
    fasta.write_text(">TX1\nAACGTCTATA\n", encoding="utf-8")
    gtf.write_text(
        'chr1\ttest\ttranscript\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "OTHER";\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="90%"):
        transcript_id_overlap(fasta, gtf, minimum=0.9)


def test_custom_upload_prepares_and_registers_transcriptome(tmp_path: Path) -> None:
    app = create_app(
        project_root=Path.cwd(),
        database_path=tmp_path / "atlas.sqlite3",
        results_dir=tmp_path / "results",
        custom_data_dir=tmp_path / "custom",
    )
    client = TestClient(app)

    response = client.post(
        "/transcriptomes/custom",
        data={
            "identifier": "uploaded-demo",
            "display_name": "Uploaded demo",
            "species": "Testus organismus",
            "assembly": "Demo1",
            "provider": "Custom",
            "release": "1",
        },
        files={
            "fasta": ("demo.fa", (DATA / "synthetic.fa").read_bytes(), "text/plain"),
            "gtf": ("demo.gtf", (DATA / "synthetic.gtf").read_bytes(), "text/plain"),
        },
        follow_redirects=False,
    )

    assert response.status_code == 303
    spec = app.state.registry.get("uploaded-demo")
    assert spec.fasta_path.exists()
    assert spec.annotation_path.exists()


def test_custom_upload_enforces_streaming_size_limit(tmp_path: Path) -> None:
    app = create_app(
        project_root=Path.cwd(),
        database_path=tmp_path / "atlas.sqlite3",
        results_dir=tmp_path / "results",
        custom_data_dir=tmp_path / "custom",
        max_upload_bytes=10,
    )
    client = TestClient(app)

    response = client.post(
        "/transcriptomes/custom",
        data={
            "identifier": "too-large",
            "display_name": "Too large",
            "species": "Testus organismus",
            "assembly": "Demo1",
            "provider": "Custom",
            "release": "1",
        },
        files={
            "fasta": ("demo.fa", b">TX1\nAACGTCTATA\n", "text/plain"),
            "gtf": ("demo.gtf", b"01234567890", "text/plain"),
        },
    )

    assert response.status_code == 413
    assert "512 MB" not in response.text
