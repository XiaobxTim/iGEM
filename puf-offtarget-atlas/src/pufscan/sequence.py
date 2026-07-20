from __future__ import annotations

from itertools import combinations, product

DNA_ALPHABET = "ACGT"


def normalize_query(query: str) -> str:
    normalized = query.strip().upper().replace("U", "T")
    if not 8 <= len(normalized) <= 12:
        raise ValueError("PUF target sequence length must be between 8 and 12 nt")
    invalid = sorted(set(normalized) - set(DNA_ALPHABET))
    if invalid:
        raise ValueError(f"PUF target sequence contains invalid characters: {', '.join(invalid)}")
    return normalized


def dna_to_rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def generate_variants(query: str, max_mismatches: int) -> dict[str, tuple[int, ...]]:
    if not 0 <= max_mismatches <= 3:
        raise ValueError("max_mismatches must be between 0 and 3 in this version")
    variants: dict[str, tuple[int, ...]] = {query: ()}
    for mismatch_count in range(1, max_mismatches + 1):
        for positions in combinations(range(len(query)), mismatch_count):
            replacements = [tuple(base for base in DNA_ALPHABET if base != query[pos]) for pos in positions]
            for bases in product(*replacements):
                variant = list(query)
                for position, base in zip(positions, bases, strict=True):
                    variant[position] = base
                variants["".join(variant)] = tuple(position + 1 for position in positions)
    return variants


def hamming_distance(query: str, observed: str) -> tuple[int, tuple[int, ...]]:
    if len(query) != len(observed):
        raise ValueError("Hamming distance requires equal-length sequences")
    positions = tuple(index for index, (left, right) in enumerate(zip(query, observed, strict=True), 1) if left != right)
    return len(positions), positions


def match_pattern(query: str, observed: str) -> str:
    if len(query) != len(observed):
        raise ValueError("Alignment requires equal-length sequences")
    bars = "".join("|" if left == right else "." for left, right in zip(query, observed, strict=True))
    return f"Query: {dna_to_rna(query)}\n       {bars}\nHit:   {dna_to_rna(observed)}"

