import pytest

from ensemble_analyzer._nmr.nuclei import (
    NUCLEUS_DATA,
    abundance,
    element_of,
    get_nucleus_spin,
    isotopes_of,
    mass_number,
    normalize_nucleus_name,
    rescale_coupling,
    shielding_to_shift,
)


@pytest.mark.parametrize("label,expected", [
    ("1H", "1H"), ("1h", "1H"), ("13-C", "13C"), ("C13", "13C"), (" 15 n ", "15N"), ("119Sn", "119Sn"),
])
def test_normalize(label, expected):
    assert normalize_nucleus_name(label) == expected


@pytest.mark.parametrize("label", ["H", "C", "12C", "foo", "999X"])
def test_normalize_rejects(label):
    with pytest.raises(ValueError):
        normalize_nucleus_name(label)


def test_data_complete():
    for key, d in NUCLEUS_DATA.items():
        assert set(d) == {"element", "spin", "gamma", "abundance", "range"}, key
        assert d["spin"] > 0 and 0 < d["abundance"] <= 1 and d["range"][0] < d["range"][1]


def test_lookups():
    assert element_of("15N") == "N" and mass_number("15N") == 15
    assert get_nucleus_spin("14N") == 1.0 and get_nucleus_spin("2H") == 1.0
    assert abundance("19F") == 1.0
    assert isotopes_of("N") == ["14N", "15N"]
    assert isotopes_of("Ag")[0] == "107Ag"


def test_shift_and_rescale():
    assert shielding_to_shift(31.8, 31.8) == 0.0
    assert shielding_to_shift(25.0, 31.8, 1.0) == pytest.approx(6.8 / (1 - 31.8e-6) + 1.0)
    # 14N -> 15N flips sign and scales by gamma ratio
    j15 = rescale_coupling(10.0, ("13C", "14N"), ("13C", "15N"))
    assert j15 == pytest.approx(10.0 * -27.116 / 19.331, rel=1e-4)
