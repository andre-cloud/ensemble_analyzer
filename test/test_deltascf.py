"""
ΔSCF: engine-neutral logic (protocol/deltascf.py), parser support and protocol
validation. Reference numbers come from ORCA 6.1.1 / NWChem 7.3.1 runs on
formaldehyde (see docs/source/user_guide/deltascf.md).
"""

from unittest.mock import MagicMock, mock_open, patch

import pytest

from ensemble_analyzer.parsers.nwchem import NWChemParser
from ensemble_analyzer.parsers.orca import OrcaParser
from ensemble_analyzer.protocol import deltascf as d
from ensemble_analyzer.protocol.deltascf import Excitation
from ensemble_analyzer.protocol.protocol import Protocol


# ---------------------------------------------------------------------------
# Excitation labels and engine translations
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text, expected", [
    ("HOMO->LUMO", (0, 0)), ("homo-1 -> lumo", (1, 0)), ("HOMO->LUMO+2", (0, 2)), ("HOMO-3->LUMO+1", (3, 1)),
])
def test_parse_excitation(text, expected):
    assert d.parse_excitation(text) == expected


@pytest.mark.parametrize("text", ["S1", "LUMO->HOMO", "HOMO+1->LUMO", ""])
def test_parse_excitation_rejects(text):
    with pytest.raises(ValueError):
        d.parse_excitation(text)


def test_orca_conf_matches_validated_lists():
    # verified against NWChem swaps on formaldehyde (same energies and <S2>)
    assert d.orca_conf(0, 0) == "0,1"        # HOMO->LUMO
    assert d.orca_conf(1, 0) == "0,1,1"      # HOMO-1->LUMO
    assert d.orca_conf(0, 1) == "0,0,1"      # HOMO->LUMO+1


def test_nwchem_swap_absolute_indices():
    assert d.nwchem_swap(0, 0, 8) == (8, 9)
    assert d.nwchem_swap(1, 0, 8) == (7, 9)
    assert d.nwchem_swap(0, 1, 8) == (8, 10)


@pytest.mark.parametrize("mult, exc, expected", [
    (1, Excitation("alpha", 0, 0), "singlet_ap"),
    (1, Excitation("beta", 2, 1), "singlet_ap"),
    (3, Excitation("beta", 0, 0), "pure"),          # T2 of formaldehyde (beta into a SOMO)
    (3, Excitation("alpha", 0, 0), "pure"),         # alpha out of a SOMO
    (3, Excitation("alpha", 1, 0), "pure"),         # second SOMO
    (3, Excitation("alpha", 2, 0), "unsupported"),  # 4 open shells
    (2, Excitation("beta", 0, 0), "pure"),          # D1 of the cation
    (2, Excitation("alpha", 0, 0), "pure"),
    (2, Excitation("alpha", 1, 0), "doublet_recoupling"),
    (2, Excitation("beta", 0, 1), "doublet_recoupling"),
])
def test_classify(mult, exc, expected):
    assert d.classify(mult, exc) == expected


# ---------------------------------------------------------------------------
# Spin purification
# ---------------------------------------------------------------------------

def test_yamaguchi_formaldehyde_s1():
    # ΔSCF mixed singlet and triplet at the optimised S1 geometry (ORCA, PBE/def2-TZVP)
    e_s = d.yamaguchi(-114.303612454, 1.006722, -114.308881638, 2.003895)
    assert e_s == pytest.approx(-114.298293, abs=1e-6)


def test_yamaguchi_reduces_to_ziegler():
    assert d.yamaguchi(-1.0, 1.0, -1.2, 2.0) == pytest.approx(2 * -1.0 - -1.2)


def test_recoupling_formaldehyde_cation():
    energies = (-113.782816951048, -113.727573897726, -113.765234807151, -113.732892073977)
    lower, upper = d.recouple_doublet(*energies)
    assert lower == pytest.approx(-113.756745, abs=1e-6)
    assert upper == pytest.approx(-113.686138, abs=1e-6)
    assert all(k > 0 for k in d.exchange_couplings(*energies))


def test_doublet_determinants_validated_case():
    # cation D0, i = HOMO-1 (pi), s = SOMO, a = LUMO (pi*): inputs used in the validation
    dets, target = d.doublet_determinants(Excitation("alpha", 1, 0))
    conf = {name: (mult, [(e.spin, d.orca_conf(e.n, e.m)) for e in exc]) for name, (mult, exc) in dets.items()}
    assert target == "M1"
    assert conf == {
        "M1": (2, [("alpha", "0,1,1")]),
        "M2": (2, [("alpha", "0,1"), ("beta", "0,1")]),
        "M3": (2, [("beta", "0,0,1")]),
        "Q": (4, []),
    }


def test_doublet_determinants_beta_target_is_m3():
    dets, target = d.doublet_determinants(Excitation("beta", 0, 1))
    assert target == "M3"
    assert dets["M3"][1] == [Excitation("beta", 0, 1)]


# ---------------------------------------------------------------------------
# Auto excitation from the reference TD-DFT output
# ---------------------------------------------------------------------------

ORCA_SINGLET = """
                         Program Version 6.1.1  -  RELEASE   -
 Number of Electrons    NEL             ....   16
 Multiplicity           Mult            ....    1
STATE  1:  E=   0.144596 au      3.935 eV    31735.2 cm**-1 <S**2> =   0.000000 Mult 1
     7a ->   8a  :     0.999843 (c= -0.99992143)

STATE  2:  E=   0.278667 au      7.583 eV    61160.3 cm**-1 <S**2> =   0.000000 Mult 1
     7a ->   9a  :     0.991597 (c= -0.99578946)

FINAL SINGLE POINT ENERGY      -114.282970
ORCA TERMINATED NORMALLY
"""

ORCA_TRIPLET = """
                         Program Version 6.1.1  -  RELEASE   -
 Number of Electrons    NEL             ....   16
 Multiplicity           Mult            ....    3
Expectation value of <S**2>     :     2.002636
STATE  1:  E=   0.111262 au      3.028 eV    24419.3 cm**-1 <S**2> =   2.003381 Mult 3
     6b ->   7b  :     0.998381 (c= -0.99918995)

STATE  2:  E=   0.145303 au      3.954 eV    31890.3 cm**-1 <S**2> =   2.002954 Mult 3
     8a ->   9a  :     0.994413 (c=  0.99720241)

FINAL SINGLE POINT ENERGY      -114.159931
ORCA TERMINATED NORMALLY
"""

NWCHEM_SINGLET = """
           Alpha electrons :     8
            Beta electrons :     8
  Root   1 singlet a2             0.144649535 a.u.                3.9361 eV
     Occ.    8  b2  ---  Virt.    9  b1       0.99992
  Root   1 triplet a2             0.100000000 a.u.                2.7000 eV
     Occ.    7  b1  ---  Virt.    9  b1       0.99000
  Root   2 singlet b2             0.278768877 a.u.                7.5857 eV
     Occ.    5  b2  ---  Virt.   10  a1       0.05907
     Occ.    8  b2  ---  Virt.   10  a1       0.99573 X
     Occ.    8  b2  ---  Virt.   10  a1       0.20000 Y
 Total times  cpu:        1.0s     wall:        1.0s
"""

NWCHEM_DOUBLET = """
           Alpha electrons :     8
            Beta electrons :     7
      <S2> =      0.7550 (Exact =     0.7500)
  Root   1         b1             0.154130913 a.u.                4.1941 eV
     Occ.    5 alpha b2  ---  Virt.    9 alpha b1      -0.05058
     Occ.    8 alpha b2  ---  Virt.    9 alpha b1      -0.20177
     Occ.    7 beta  b1  ---  Virt.    8 beta  b2      -0.97739
  Root   2         a1             0.198830917 a.u.                5.4105 eV
     Occ.    6 beta  a1  ---  Virt.    8 beta  b2       0.99777
 Total times  cpu:        1.0s     wall:        1.0s
"""


def _parser(cls, text):
    with patch("builtins.open", mock_open(read_data=text)):
        return cls("dummy.out", MagicMock())


def _auto(parser, mult, root, base):
    n_alpha, n_beta = parser.parse_n_electrons()
    return d.resolve_auto(parser.parse_excited_states(mult), root, n_alpha, n_beta, base)


def test_orca_singlet_auto():
    p = _parser(OrcaParser, ORCA_SINGLET)
    assert p.parse_n_electrons() == (8, 8)
    assert _auto(p, 1, 1, 0) == Excitation("alpha", 0, 0, weight=0.999843, source="auto")
    assert _auto(p, 1, 2, 0).label == "HOMO->LUMO+1 (alpha)"


def test_orca_triplet_auto_and_s2():
    p = _parser(OrcaParser, ORCA_TRIPLET)
    assert p.parse_n_electrons() == (9, 7)
    assert p.parse_s2() == pytest.approx(2.002636)
    assert _auto(p, 3, 1, 0).label == "HOMO->LUMO (beta)"     # T2 = beta 6b -> 7b
    assert _auto(p, 3, 2, 0).label == "HOMO->LUMO (alpha)"    # T3 = alpha 8a -> 9a
    assert [s.energy_ev for s in p.parse_excited_states(3)] == [3.028, 3.954]


def test_nwchem_singlet_auto_skips_triplets_and_y():
    p = _parser(NWChemParser, NWCHEM_SINGLET)
    states = p.parse_excited_states(1)
    assert [s.root for s in states] == [1, 2]
    assert states[0].energy_ev == pytest.approx(3.9361)
    assert _auto(p, 1, 1, 1).label == "HOMO->LUMO (alpha)"
    assert _auto(p, 1, 2, 1).label == "HOMO->LUMO+1 (alpha)"
    assert all(t[3] < 1.0 for t in states[1].transitions)       # Y line skipped


def test_nwchem_doublet_auto():
    p = _parser(NWChemParser, NWCHEM_DOUBLET)
    assert p.parse_n_electrons() == (8, 7)
    assert p.parse_s2() == pytest.approx(0.7550)
    first = _auto(p, 2, 1, 1)
    assert first.label == "HOMO->LUMO (beta)"
    assert first.weight == pytest.approx(0.97739 ** 2)
    assert _auto(p, 2, 2, 1).label == "HOMO-1->LUMO (beta)"


def test_resolve_auto_missing_root():
    p = _parser(OrcaParser, ORCA_SINGLET)
    with pytest.raises(ValueError, match="root 5 not found"):
        _auto(p, 1, 5, 0)


def test_orca_status_strings():
    failed = ORCA_SINGLET.replace("ORCA TERMINATED NORMALLY", "SCF NOT CONVERGED AFTER   2 CYCLES")
    assert _parser(OrcaParser, failed).deltascf_status() == "SCF not converged"
    assert _parser(OrcaParser, ORCA_SINGLET).deltascf_status() is None
    assert _parser(OrcaParser, ORCA_SINGLET).deltascf_status(opt=True) == "geometry optimisation not converged"
    bug = ORCA_SINGLET + "Error (TGeneralVectorSet): Constructor called with NVecs<=0. This is a bug\n"
    assert _parser(OrcaParser, bug).deltascf_status() == "ORCA internal error"


def test_nwchem_status_strings():
    assert _parser(NWChemParser, NWCHEM_SINGLET).deltascf_status() is None
    failed = NWCHEM_SINGLET + "      Calculation failed to converge\n"
    assert _parser(NWChemParser, failed).deltascf_status() == "SCF not converged"


# ---------------------------------------------------------------------------
# Settings and protocol validation
# ---------------------------------------------------------------------------

def _p(number, **kw):
    base = dict(functional="PBE", basis="def2-SVP", calculator="orca")
    base.update(kw)
    return Protocol(number=number, **base)


def test_settings_defaults_and_normalisation():
    s = _p(1, deltascf={"excitation": "homo-1 -> lumo+1", "spin": "Beta"}).deltascf
    assert s["excitation"] == "HOMO-1->LUMO+1" and s["spin"] == "beta"
    assert s["spin_purify"] is True and s["min_gap_ev"] == 0.1 and s["cation_guess"] is True
    assert _p(2, deltascf={}).deltascf["excitation"] == "auto"
    assert _p(3, deltascf=True).deltascf["root"] == 1
    assert _p(4).deltascf is None


@pytest.mark.parametrize("kw, match", [
    (dict(deltascf={"foo": 1}), "unknown deltascf option"),
    (dict(deltascf={"excitation": "S1"}), "not understood"),
    (dict(deltascf={"spin": "up"}), "spin must be"),
    (dict(deltascf={"alphaconf": "0,1"}, calculator="nwchem"), "ORCA-only"),
    (dict(deltascf={"alphaconf": "zero"}), "not a list of occupations"),
    (dict(deltascf={}, nroots=5), "cannot also request TD-DFT"),
    (dict(deltascf={}, add_input="! DeltaSCF UKS"), "remove"),
    (dict(deltascf={}, calculator="tblite", functional="GFN2-xTB"), "available for"),
])
def test_settings_rejected(kw, match):
    with pytest.raises(ValueError, match=match):
        _p(5, **kw)


def test_orca_older_than_6_rejected(monkeypatch):
    monkeypatch.setenv("ORCAVERSION", "5.0.4")
    with pytest.raises(ValueError, match="ORCA >= 6.0"):
        _p(1, deltascf={})


def test_ladder_selection(monkeypatch):
    monkeypatch.setenv("ORCAVERSION", "6.1.1")
    assert d.ladder("orca", {})[0].keywords == "GMF"
    assert all(a.guess == "reference" for a in d.ladder("orca", {"cation_guess": False}))
    assert len(d.ladder("orca", {"max_attempts": 2})) == 2
    monkeypatch.setenv("ORCAVERSION", "6.0.1")
    assert all("GMF" not in a.keywords and "FreezeAndRelease" not in a.keywords for a in d.ladder("orca", {}))
    assert d.ladder("nwchem", {})[0].dft == {"max_ovl": None}


def test_check_protocols_reference():
    ref = _p(0, nroots=5, tda=True)
    same = _p(1, read_orbitals="0", deltascf={})
    auto = _p(2, deltascf={"excitation": "HOMO->LUMO"})
    log = MagicMock()
    d.check_protocols([ref, same, auto], log)
    assert same.deltascf["reference_step"] == 0
    assert auto.deltascf["reference_step"] is None

    other_level = _p(3, functional="B3LYP", read_orbitals="0", deltascf={})
    d.check_protocols([ref, other_level], log)
    assert other_level.deltascf["reference_step"] is None
    log.warning.assert_called()


@pytest.mark.parametrize("steps, match", [
    (lambda: [_p(0, nroots=5), _p(1, read_orbitals="0", deltascf={}), _p(2, read_orbitals="1", deltascf={})],
     "points to the deltascf step"),
    (lambda: [_p(0), _p(1, read_orbitals="0", deltascf={})], "must be a TD-DFT step"),
    (lambda: [_p(1, read_orbitals="7", deltascf={})], "is not a protocol step"),
])
def test_check_protocols_rejects(steps, match):
    with pytest.raises(ValueError, match=match):
        d.check_protocols(steps(), MagicMock())


# ---------------------------------------------------------------------------
# "inherit": ΔSCF opt -> ΔSCF single points (functional benchmark) without TD-DFT
# ---------------------------------------------------------------------------

def test_inherit_chain_and_explicit_source():
    a = _p(0, opt=True, deltascf={})
    b = _p(1, functional="PBE0", deltascf={"excitation": "inherit"})
    c = _p(2, functional="CAM-B3LYP", deltascf={"excitation": "Inherit"})
    e = _p(3, functional="B3LYP", deltascf={"excitation": "inherit", "inherit_from": 0})
    d.check_protocols([a, b, c, e], MagicMock())
    assert (b.deltascf["inherit_step"], c.deltascf["inherit_step"], e.deltascf["inherit_step"]) == (0, 1, 0)
    assert b.deltascf["inherit_gap_tol_ev"] == 1.0


@pytest.mark.parametrize("steps, match", [
    (lambda: [_p(1, deltascf={"excitation": "inherit"})], "needs an earlier deltascf step"),
    (lambda: [_p(0, mult=3, deltascf={}), _p(1, deltascf={"excitation": "inherit"})], "same multiplicity"),
    (lambda: [_p(0, deltascf={}), _p(1, deltascf={"excitation": "inherit", "inherit_from": 4})], "inherit_from=4"),
])
def test_inherit_source_errors(steps, match):
    with pytest.raises(ValueError, match=match):
        d.check_protocols(steps(), MagicMock())


@pytest.mark.parametrize("kw, match", [
    ({"excitation": "inherit", "alphaconf": "0,1"}, "cannot be combined"),
    ({"inherit_from": 0}, "requires excitation 'inherit'"),
    ({"excitation": "inherit", "inherit_gap_tol_ev": -1}, "inherit_gap_tol_ev"),
])
def test_inherit_settings_rejected(kw, match):
    with pytest.raises(ValueError, match=match):
        _p(5, deltascf=kw)
