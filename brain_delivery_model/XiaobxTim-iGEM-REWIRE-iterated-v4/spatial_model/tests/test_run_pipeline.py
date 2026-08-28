from __future__ import annotations

import pytest

from spatial_model.run_pipeline import parse_doses, scenario_name


def test_parse_doses_has_stable_three_scenario_contract():
    assert parse_doses("0.3,1.0,3.0") == (0.3, 1.0, 3.0)
    assert [scenario_name(value) for value in parse_doses("0.3,1.0,3.0")] == [
        "dose_0p3",
        "dose_1p0",
        "dose_3p0",
    ]


@pytest.mark.parametrize("value", ["", "1.0,1.0", "-1,1", "nan,1", "0.3,1.0"])
def test_parse_doses_rejects_unsafe_comparison_sets(value):
    with pytest.raises(ValueError):
        parse_doses(value)
