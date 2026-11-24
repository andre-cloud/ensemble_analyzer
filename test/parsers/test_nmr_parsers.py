"""NMR parsing on real outputs (ORCA 6.1, Gaussian 16 A.03, NWChem 7), trimmed to the NMR sections."""

from pathlib import Path

import pytest

from ensemble_analyzer._nmr.parsing import parse_nmr_output
from ensemble_analyzer.protocol.nmr import normalize_settings

DATA = Path(__file__).parent / "data" / "nmr"
FLUOROETHANE = ["C", "C", "F", "H", "H", "H", "H", "H"]
ACETONITRILE = ["C", "C", "N", "H", "H", "H"]
TMS = ["C", "Si", "C", "C", "C"] + ["H"] * 12


def _settings(engine, nuclei, atoms=None, shielding=True, coupling=True):
    raw = {"shielding": shielding, "coupling": coupling, "nuclei": nuclei, "atoms": atoms,
           "references": [{"name": "x", "smiles": "C", "nuclei": nuclei}] if shielding else []}
    return normalize_settings(raw, 1, engine)


def _parse(engine, key, symbols, *nuclei, **kw):
    text = (DATA / f"{engine}_{key}.out").read_text()
    return parse_nmr_output(engine, text, _settings(engine, list(nuclei), **kw), symbols)


# --- ORCA -------------------------------------------------------------------

def test_orca_full():
    r = _parse("orca", "A", FLUOROETHANE, "1H", "13C", "19F")
    assert r.shielding[0] == pytest.approx(173.550) and r.shielding[2] == pytest.approx(422.852)
    assert r.aniso[1] == pytest.approx(82.465)
    assert len(r.coupling) == 28
    assert r.j(0, 1) == pytest.approx(96.316, abs=1e-3)
    assert r.j(1, 2) == pytest.approx(-274.944, abs=1e-3)
    assert r.j(3, 4) == pytest.approx(-11.837, abs=1e-3)


def test_orca_subset_with_isotopes():
    r = _parse("orca", "E", FLUOROETHANE, "1H", "13C", "19F", atoms=[1, 2, 6, 7])
    assert set(r.shielding) == {1, 2, 6, 7}
    assert r.j(1, 2) == pytest.approx(-274.570, abs=1e-3)
    assert r.j(1, 6) == pytest.approx(196.712, abs=1e-3)


def test_orca_pairs_without_isotope_are_dropped():
    # atom list without ist: ORCA prints "1C gnA = 0.000" and J = 0 for every heteronuclear pair
    r = _parse("orca", "B", FLUOROETHANE, "1H", "13C", "19F", atoms=[1, 2, 6, 7])
    assert set(r.coupling) == {(6, 7)}
    assert r.j(6, 7) == pytest.approx(-7.224, abs=1e-3)


def test_orca_isotope_rescaling():
    r15 = _parse("orca", "D", ACETONITRILE, "1H", "13C", "15N")
    assert r15.j(1, 2) == pytest.approx(-51.486, abs=1e-3)
    r14 = _parse("orca", "D", ACETONITRILE, "1H", "13C", "14N")
    assert r14.j(1, 2) == pytest.approx(-51.486 * 19.331 / -27.116, rel=1e-4)


def test_orca_shielding_only_keeps_requested_nuclei():
    r = _parse("orca", "C", TMS, "1H", coupling=False)
    assert set(r.shielding) == set(range(5, 17)) and not r.coupling


# --- Gaussian ---------------------------------------------------------------

def test_gaussian_full_from_reduced_coupling():
    r = _parse("gaussian", "A", FLUOROETHANE, "1H", "13C", "19F")
    assert r.shielding[0] == pytest.approx(173.5776) and r.aniso[2] == pytest.approx(56.8603)
    assert r.j(0, 1) == pytest.approx(95.7592, abs=1e-3)
    assert r.j(1, 2) == pytest.approx(-276.973, abs=1e-2)


def test_gaussian_15n_matches_printed_j():
    # Gaussian uses gamma(15N) from mu = -0.28319 (0.04 % from the IUPAC value in nuclei.json)
    r = _parse("gaussian", "D", ACETONITRILE, "1H", "13C", "15N")
    assert r.j(1, 2) == pytest.approx(-50.7982, rel=1e-3)
    assert r.j(0, 1) == pytest.approx(154.307, abs=1e-2)


def test_gaussian_subset_filters_atoms():
    r = _parse("gaussian", "B", FLUOROETHANE, "1H", "13C", "19F", atoms=[1, 2, 6, 7])
    assert set(r.shielding) == {1, 2, 6, 7}
    assert set(r.coupling) <= {(1, 2), (1, 6), (1, 7), (2, 6), (2, 7), (6, 7)}
    assert r.j(1, 2) == pytest.approx(-276.973, abs=1e-2)


# --- NWChem -----------------------------------------------------------------

def test_nwchem_full():
    r = _parse("nwchem", "E", FLUOROETHANE, "1H", "13C", "19F")
    assert r.shielding[0] == pytest.approx(173.5575)
    assert r.j(0, 1) == pytest.approx(95.9924, abs=1e-3)
    assert r.j(0, 2) == pytest.approx(55.0519, abs=1e-3)


def test_nwchem_default_14n_rescaled_to_15n():
    r = _parse("nwchem", "F", ACETONITRILE, "1H", "13C", "15N")
    assert r.j(0, 2) == pytest.approx(-2.9166 * -27.116 / 19.331, rel=1e-4)


def test_nwchem_subset_shielding_is_renumbered():
    r = _parse("nwchem", "B", FLUOROETHANE, "1H", "13C", "19F", atoms=[1, 2, 6, 7], coupling=False)
    assert r.shielding == pytest.approx({1: 116.8568, 2: 423.0989, 6: 27.4001, 7: 27.4001})


def test_nwchem_tms_shielding():
    r = _parse("nwchem", "C", TMS, "1H", "13C", coupling=False)
    assert r.shielding[0] == pytest.approx(191.0721) and 1 not in r.shielding
