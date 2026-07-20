from pathlib import Path

import pytest

from pufscan.structure import RNAplfoldAnalyzer, parse_lunp


def test_parse_lunp_extracts_per_base_and_motif_probability(tmp_path: Path) -> None:
    lunp = tmp_path / "plfold_lunp"
    lunp.write_text(
        "#i$\tl=1\t2\t3\n1\t0.8\t0.5\t0.2\n2\t0.6\t0.4\t0.1\n3\t0.4\t0.3\t0.05\n",
        encoding="utf-8",
    )
    table = parse_lunp(lunp)
    assert table[1][0] == pytest.approx(0.8)
    assert table[2][2] == pytest.approx(0.1)


def test_missing_rnaplfold_returns_na_without_random_values(tmp_path: Path) -> None:
    analyzer = RNAplfoldAnalyzer(cache_dir=tmp_path, executable="definitely-not-installed-RNAplfold")
    result = analyzer.analyze("AACGTCTATA", 0, 10)
    assert result.motif_mean_unpaired_probability is None
    assert result.accessibility_score is None
    assert "not available" in result.warnings[0]

