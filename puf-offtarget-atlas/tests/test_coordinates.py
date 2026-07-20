from pathlib import Path

from pufscan.coordinates import TranscriptCoordinateIndex

DATA = Path(__file__).parent / "data"


def test_positive_strand_hit_spans_exon_junction() -> None:
    index = TranscriptCoordinateIndex.from_gtf(DATA / "synthetic.gtf")
    mapped = index.map_hit("ENST000001.1", 10, 20, splice_proximity_nt=20)
    assert mapped.genomic_blocks_1based == (("chr1", 111, 115), ("chr1", 201, 205))
    assert mapped.junction_spanning is True
    assert mapped.transcript_region == "splice_junction_spanning"
    assert mapped.distance_to_nearest_junction == 0


def test_negative_strand_coordinates_follow_transcript_orientation() -> None:
    index = TranscriptCoordinateIndex.from_gtf(DATA / "synthetic.gtf")
    mapped = index.map_hit("ENST000002.1", 5, 15, splice_proximity_nt=2)
    assert mapped.genomic_blocks_1based == (("chr2", 401, 410),)
    assert mapped.strand == "-"
    assert mapped.exon_numbers == ("1",)


def test_generic_utr_is_split_using_cds_bounds() -> None:
    index = TranscriptCoordinateIndex.from_gtf(DATA / "synthetic.gtf")
    five_prime = index.map_hit("ENST000004.1", 0, 10)
    three_prime = index.map_hit("ENST000004.1", 30, 40)
    assert five_prime.transcript_region == "5UTR"
    assert three_prime.transcript_region == "3UTR"


def test_splice_proximal_hit_is_flagged_without_crossing() -> None:
    index = TranscriptCoordinateIndex.from_gtf(DATA / "synthetic.gtf")
    mapped = index.map_hit("ENST000001.1", 0, 8, splice_proximity_nt=8)
    assert mapped.junction_spanning is False
    assert mapped.transcript_region == "splice_proximal"
    assert mapped.distance_to_nearest_junction == 7

