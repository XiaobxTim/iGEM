from pathlib import Path

from pufscan.search import SearchEngine, TranscriptRecord, iter_transcripts

DATA = Path(__file__).parent / "data"


def test_aho_search_finds_exact_and_one_mismatch_hits() -> None:
    result = SearchEngine("AACGUCUAUA", max_mismatches=1).search(iter_transcripts(DATA / "synthetic.fa"))
    assert result.statistics.transcripts_scanned == 5
    assert result.statistics.total_nucleotides_scanned == 165
    assert len(result.hits) == 6
    assert sum(hit.mismatch_count == 0 for hit in result.hits) == 5
    mismatch = next(hit for hit in result.hits if hit.transcript_id == "ENST000005.1")
    assert mismatch.mismatch_positions_1based == (5,)
    assert mismatch.matched_sequence_rna == "AACGACUAUA"


def test_search_respects_mismatch_threshold() -> None:
    result = SearchEngine("AACGUCUAUA", max_mismatches=0).search(iter_transcripts(DATA / "synthetic.fa"))
    assert all(hit.mismatch_count == 0 for hit in result.hits)
    assert all(hit.transcript_id != "ENST000005.1" for hit in result.hits)


def test_reverse_complement_is_opt_in_and_labeled() -> None:
    record = TranscriptRecord(
        transcript_id="TX.1",
        transcript_id_without_version="TX",
        gene_id="GENE.1",
        gene_id_without_version="GENE",
        gene_name="TEST",
        transcript_type="synthetic",
        sequence="CCCTATAGACGTTCCC",
    )
    forward = SearchEngine("AACGUCUAUA", max_mismatches=0).search([record])
    reverse = SearchEngine("AACGUCUAUA", max_mismatches=0, search_reverse_complement=True).search([record])
    assert forward.hits == ()
    assert reverse.hits[0].search_orientation == "reverse_complement"


def test_twelve_nt_query_is_supported() -> None:
    record = TranscriptRecord("TX", "TX", None, None, None, None, "GGAACGTCTATACC")
    result = SearchEngine("AACGUCUAUACC", max_mismatches=0).search([record])
    assert len(result.hits) == 1

