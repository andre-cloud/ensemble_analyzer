"""
TD-DFT input generation and energy selection.

- ORCA: ``iroot 0`` is added to the %tddft block of spectra steps, so that
  FINAL SINGLE POINT ENERGY is the ground state.
- ``_wants_excited_state_energy``: when the parsed energy must stay the one of
  the followed excited root.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from ensemble_analyzer._parser_parameter import _wants_excited_state_energy
from ensemble_analyzer.calculators.orca import OrcaCalc


@pytest.fixture
def tddft_protocol(mock_protocol):
    mock_protocol.charge = 0
    mock_protocol.mult = 1
    mock_protocol.solvent = None
    mock_protocol.add_input = ""
    mock_protocol.read_orbitals = ""
    mock_protocol.constraints = []
    mock_protocol.ts = False
    mock_protocol.opt = False
    mock_protocol.freq = False
    mock_protocol.nroots = 10
    mock_protocol.tda = True
    return mock_protocol


@pytest.fixture
def conf(mock_conformer):
    mock_conformer.folder = "conf_1"
    mock_conformer.number = 1
    return mock_conformer


def _orca_tddft_block(protocol, conformer, method="single_point") -> str:
    with patch("enan_calculators._orca.shutil.which", return_value="orca"), \
         patch("ensemble_analyzer.calculators.orca.orca_profile"), \
         patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 6):
        ase_calc, _ = getattr(OrcaCalc(protocol, 4, conformer), method)()
    blocks = ase_calc.parameters["orcablocks"]
    # the generated block (the user's own %tddft block from add_input may precede it)
    return next(line for line in blocks.splitlines() if line.startswith("%tddft nroots"))


class TestOrcaIroot:

    def test_single_point_follows_ground_state(self, tddft_protocol, conf):
        assert _orca_tddft_block(tddft_protocol, conf) == "%tddft nroots 10 tda true iroot 0 end"

    def test_optimisation_keeps_root(self, tddft_protocol, conf):
        tddft_protocol.opt = True
        block = _orca_tddft_block(tddft_protocol, conf, "optimisation")
        assert "iroot" not in block

    def test_opt_in_add_input_keeps_root(self, tddft_protocol, conf):
        tddft_protocol.add_input = "! Opt"
        assert "iroot" not in _orca_tddft_block(tddft_protocol, conf)

    def test_user_iroot_not_overridden(self, tddft_protocol, conf):
        tddft_protocol.add_input = "%tddft iroot 2 end"
        assert "iroot" not in _orca_tddft_block(tddft_protocol, conf)

    def test_full_tddft_flag(self, tddft_protocol, conf):
        tddft_protocol.tda = False
        assert _orca_tddft_block(tddft_protocol, conf) == "%tddft nroots 10 tda false iroot 0 end"


class TestExcitedStateEnergyFlag:

    @staticmethod
    def _p(**kw):
        base = dict(nroots=10, opt=False, add_input="")
        base.update(kw)
        return SimpleNamespace(**base)

    def test_spectra_step_uses_ground_state(self):
        assert _wants_excited_state_energy(self._p()) is False

    def test_no_tddft(self):
        assert _wants_excited_state_energy(self._p(nroots=None, opt=True)) is False

    def test_excited_state_optimisation(self):
        assert _wants_excited_state_energy(self._p(opt=True)) is True
        assert _wants_excited_state_energy(self._p(add_input="! Opt")) is True

    def test_explicit_iroot(self):
        assert _wants_excited_state_energy(self._p(add_input="%tddft IRoot 2 end")) is True
        assert _wants_excited_state_energy(self._p(add_input="%tddft iroot 0 end")) is False

    def test_mock_protocol_is_safe(self):
        assert _wants_excited_state_energy(MagicMock()) is False
