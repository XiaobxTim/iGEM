from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from Bio.Seq import Seq

if TYPE_CHECKING:
    from pufscan.coordinates import TranscriptCoordinateIndex

Editor = Literal["APOBEC_C2U", "ADAR_A2I"]


def parse_editing_window(raw: str) -> tuple[int, int]:
    try:
        start_raw, end_raw = raw.split(":", 1)
        start, end = int(start_raw), int(end_raw)
    except (ValueError, TypeError) as error:
        raise ValueError("Editing window must use START:END integer offsets") from error
    if start > end:
        raise ValueError("Editing window start must not exceed end")
    return start, end


@dataclass(frozen=True)
class CodingEditResult:
    reference_codon: str | None
    edited_codon: str | None
    reference_amino_acid: str | None
    edited_amino_acid: str | None
    protein_position: int | None
    coding_consequence: str


@dataclass(frozen=True)
class PotentialEditingEvent:
    transcript_id: str
    editor: Editor
    event_label: str
    reference_base: str
    simulated_base: str
    position_relative_to_motif_start: int
    position_relative_to_motif_end: int
    transcript_position_1based: int
    chromosome: str
    genomic_position_1based: int
    transcript_region: str
    reference_codon: str | None
    edited_codon: str | None
    reference_amino_acid: str | None
    edited_amino_acid: str | None
    protein_position: int | None
    coding_consequence: str
    consequence_summary: str
    consequence_evidence_level: str


def classify_coding_edit(
    transcript_sequence: str,
    cds_positions0: tuple[int, ...],
    edit_position0: int,
    editor: Editor,
) -> CodingEditResult:
    replacement = {"APOBEC_C2U": ("C", "T"), "ADAR_A2I": ("A", "G")}[editor]
    reference_base, edited_base = replacement
    sequence = transcript_sequence.upper().replace("U", "T")
    if edit_position0 < 0 or edit_position0 >= len(sequence):
        raise ValueError("Edit position lies outside transcript sequence")
    if sequence[edit_position0] != reference_base:
        raise ValueError(f"Base at edit position is not editable by {editor}")
    try:
        cds_offset = cds_positions0.index(edit_position0)
    except ValueError:
        return CodingEditResult(None, None, None, None, None, "unknown")
    codon_start = cds_offset - cds_offset % 3
    codon_positions = cds_positions0[codon_start : codon_start + 3]
    if len(codon_positions) != 3:
        return CodingEditResult(None, None, None, None, None, "unknown")
    reference_codon = "".join(sequence[position] for position in codon_positions)
    edited = list(reference_codon)
    edited[cds_offset % 3] = edited_base
    edited_codon = "".join(edited)
    reference_amino_acid = str(Seq(reference_codon).translate())  # type: ignore[no-untyped-call]
    edited_amino_acid = str(Seq(edited_codon).translate())  # type: ignore[no-untyped-call]
    protein_position = codon_start // 3 + 1
    if reference_amino_acid != "*" and edited_amino_acid == "*":
        consequence = "stop_gained"
    elif reference_amino_acid == "*" and edited_amino_acid != "*":
        consequence = "stop_lost"
    elif protein_position == 1 and reference_amino_acid == "M" and edited_amino_acid != "M":
        consequence = "start_lost"
    elif reference_amino_acid == edited_amino_acid:
        consequence = "synonymous"
    else:
        consequence = "missense"
    return CodingEditResult(
        reference_codon,
        edited_codon,
        reference_amino_acid,
        edited_amino_acid,
        protein_position,
        consequence,
    )


def find_potential_editing_events(
    transcript_id: str,
    transcript_sequence: str,
    motif_start0: int,
    motif_end0: int,
    editor: Editor,
    editing_window: tuple[int, int],
    coordinate_index: TranscriptCoordinateIndex,
    splice_proximity_nt: int = 20,
) -> tuple[PotentialEditingEvent, ...]:
    from pufscan.coordinates import TranscriptCoordinateIndex

    if not isinstance(coordinate_index, TranscriptCoordinateIndex):
        raise TypeError("coordinate_index must be a TranscriptCoordinateIndex")
    reference_base, simulated_base = {
        "APOBEC_C2U": ("C", "T"),
        "ADAR_A2I": ("A", "G"),
    }[editor]
    sequence = transcript_sequence.upper().replace("U", "T")
    window_start = max(0, motif_start0 + editing_window[0])
    window_end = min(len(sequence), motif_start0 + editing_window[1] + 1)
    annotation = coordinate_index.get(transcript_id)
    cds_positions = tuple(
        position
        for start0, end0 in sorted(annotation.feature_intervals.get("CDS", []))
        for position in range(start0, end0)
    )
    events: list[PotentialEditingEvent] = []
    for position0 in range(window_start, window_end):
        if sequence[position0] != reference_base:
            continue
        mapped = coordinate_index.map_hit(
            transcript_id, position0, position0 + 1, splice_proximity_nt=splice_proximity_nt
        )
        if position0 in cds_positions:
            coding = classify_coding_edit(sequence, cds_positions, position0, editor)
            summary = f"Potential {coding.coding_consequence} coding change. Experimental validation is required."
            level = "Level 1"
        else:
            coding = CodingEditResult(None, None, None, None, None, "unknown")
            if mapped.transcript_region == "5UTR":
                summary = "Potential effect on translation initiation or upstream regulatory elements. Experimental validation is required."
            elif mapped.transcript_region == "3UTR":
                summary = "Potential effect on RNA stability, localization, miRNA binding or RBP binding. No specific regulatory disruption is claimed without external evidence."
            elif mapped.transcript_region in {"splice_proximal", "splice_junction_spanning"}:
                summary = "Potential effect on RNA splicing or exon definition."
            else:
                summary = "Potential alteration of RNA structure or regulatory function."
            level = "Level 1"
        genomic_position = mapped.genomic_blocks_1based[0][1]
        events.append(
            PotentialEditingEvent(
                transcript_id=transcript_id,
                editor=editor,
                event_label="potential editable base",
                reference_base=reference_base,
                simulated_base=simulated_base,
                position_relative_to_motif_start=position0 - motif_start0,
                position_relative_to_motif_end=position0 - (motif_end0 - 1),
                transcript_position_1based=position0 + 1,
                chromosome=mapped.chromosome,
                genomic_position_1based=genomic_position,
                transcript_region=mapped.transcript_region,
                reference_codon=coding.reference_codon,
                edited_codon=coding.edited_codon,
                reference_amino_acid=coding.reference_amino_acid,
                edited_amino_acid=coding.edited_amino_acid,
                protein_position=coding.protein_position,
                coding_consequence=coding.coding_consequence,
                consequence_summary=summary,
                consequence_evidence_level=level,
            )
        )
    return tuple(events)
