import pytest

from pufscan.sequence import (
    generate_variants,
    hamming_distance,
    match_pattern,
    normalize_query,
    reverse_complement,
)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("aacgucua", "AACGTCTA"), ("AACGTCTA", "AACGTCTA"), (" aacgucua\n", "AACGTCTA")],
)
def test_normalize_query_accepts_rna_and_dna(raw: str, expected: str) -> None:
    assert normalize_query(raw) == expected


@pytest.mark.parametrize("query", ["ACGU", "ACGUACGUACGUA", "AACGXCUA"])
def test_normalize_query_rejects_invalid_input(query: str) -> None:
    with pytest.raises(ValueError):
        normalize_query(query)


def test_generate_variants_tracks_mismatch_positions() -> None:
    variants = generate_variants("AACGTCTA", 1)
    assert len(variants) == 25
    assert variants["AACGTCTA"] == ()
    assert variants["CACGTCTA"] == (1,)


def test_generate_variants_rejects_unsupported_mismatch_count() -> None:
    with pytest.raises(ValueError, match="0 and 3"):
        generate_variants("AACGTCTA", 4)


def test_hamming_and_alignment_are_explainable() -> None:
    query = "AACGTCTA"
    hit = "AACGACTA"
    assert hamming_distance(query, hit) == (1, (5,))
    assert match_pattern(query, hit) == "Query: AACGUCUA\n       ||||.|||\nHit:   AACGACUA"


def test_reverse_complement_uses_dna_alphabet() -> None:
    assert reverse_complement("AACGTCTA") == "TAGACGTT"

