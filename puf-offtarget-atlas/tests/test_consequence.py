from pufscan.consequence import classify_coding_edit, parse_editing_window


def test_parse_editing_window_is_inclusive() -> None:
    assert parse_editing_window("-15:10") == (-15, 10)


def test_apobec_c_to_u_can_be_synonymous() -> None:
    result = classify_coding_edit("ATGGCC", tuple(range(6)), 5, "APOBEC_C2U")
    assert result.reference_codon == "GCC"
    assert result.edited_codon == "GCT"
    assert result.coding_consequence == "synonymous"


def test_apobec_c_to_u_can_be_missense() -> None:
    result = classify_coding_edit("ATGTCC", tuple(range(6)), 4, "APOBEC_C2U")
    assert result.edited_codon == "TTC"
    assert result.coding_consequence == "missense"


def test_apobec_c_to_u_can_gain_stop() -> None:
    result = classify_coding_edit("ATGCAA", tuple(range(6)), 3, "APOBEC_C2U")
    assert result.edited_codon == "TAA"
    assert result.coding_consequence == "stop_gained"


def test_adar_a_to_i_is_read_as_g_and_can_lose_start() -> None:
    result = classify_coding_edit("ATGGCC", tuple(range(6)), 0, "ADAR_A2I")
    assert result.edited_codon == "GTG"
    assert result.coding_consequence == "start_lost"


def test_non_editable_base_is_rejected() -> None:
    try:
        classify_coding_edit("ATGGCC", tuple(range(6)), 1, "APOBEC_C2U")
    except ValueError as error:
        assert "not editable" in str(error)
    else:
        raise AssertionError("Expected a non-editable base error")

