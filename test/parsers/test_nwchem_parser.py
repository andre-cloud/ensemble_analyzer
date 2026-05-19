import pytest
import numpy as np
from unittest.mock import MagicMock, mock_open, patch
from ensemble_analyzer._parsers._nwchem import NWChemParser


NWCHEM_OPT = """
Output coordinates
  O    0.00000000    0.00000000    0.00000000
  H    0.00000000    0.00000000    1.80000000

Total DFT energy = -76.12345678

Rotational constants (GHz): A = 1.23456 B = 0.98765 C = 0.87654

Dipole Moment (Debye): X = 0.0000 Y = 0.0000 Z = 1.8500

Optimization converged

Total times  cpu:      10.0s     wall:      12.0s
"""

NWCHEM_FREQ = """
 ----------------------------------------------------------------------------
 Normal Eigenvalue ||           Projected Infra Red Intensities
  Mode   [cm**-1]  || [atomic units] [(debye/angs)**2] [(KM/mol)] [arbitrary]
 ------ ---------- || -------------- ----------------- ---------- -----------
     1      100.00 ||    0.003779           0.087         3.684       0.750
     2      200.00 ||    0.119336           2.753       116.334      23.682
     3      300.00 ||    0.119331           2.753       116.330      23.681
 ----------------------------------------------------------------------------


Total times  cpu:      10.0s     wall:      12.0s
"""

NWCHEM_TDDFT = """
Excitation energies
       1:     5.0000    0.1234
       2:     6.0000    0.5678



Total times  cpu:      10.0s     wall:      12.0s
"""

NWCHEM_CRASHED = """
Some partial output
"""

NWCHEM_NO_GEOM = """
Total DFT energy = -76.12345678

Total times  cpu:      10.0s     wall:      12.0s
"""


class TestNWChemParser:

    @pytest.fixture
    def parser(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            p.fl = ""
            return p

    def test_normal_termination(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            assert p.normal_termination() is True

    def test_crashed_calculation(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_CRASHED)):
            p = NWChemParser("dummy.log", mock_logger)
            assert p.correct_exiting is False
            assert p.normal_termination() is False
            p.log.warning.assert_called()

    def test_parse_geom(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            geom = p.parse_geom()
            assert geom.shape == (2, 3)
            assert np.allclose(geom[0], [0.0, 0.0, 0.0])
            assert np.allclose(geom[1], [0.0, 0.0, 1.8])

    def test_parse_geom_not_found(self, parser):
        parser.fl = "no geometry here"
        geom = parser.parse_geom()
        assert len(geom) == 0

    def test_parse_energy(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            assert p.parse_energy() == -76.12345678

    def test_parse_energy_not_found(self, parser):
        parser.fl = "no energy here"
        assert parser.parse_energy() == 0.0

    def test_opt_done(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            assert p.opt_done() is True

    def test_opt_not_done(self, parser):
        parser.fl = "no optimization"
        assert parser.opt_done() is False

    def test_parse_B_m(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT)):
            p = NWChemParser("dummy.log", mock_logger)
            B, M = p.parse_B_m()
            assert np.allclose(B, [1.23456 / 29.979, 0.98765 / 29.979, 0.87654 / 29.979], atol=1e-4)
            assert np.allclose(M, [0.0, 0.0, 1.85])

    def test_parse_B_m_not_found(self, parser):
        parser.fl = "no B or M"
        B, M = parser.parse_B_m()
        assert np.allclose(B, [1, 0, 0])
        assert np.allclose(M, [1, 0, 0])

    def test_parse_freq(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT + "\n" + NWCHEM_FREQ)):
            p = NWChemParser("dummy.log", mock_logger)
            freq, ir, vcd = p.parse_freq()
            assert len(freq) == 3
            assert np.allclose(freq, [100.0, 200.0, 300.0])
            assert np.allclose(ir[:, 0], [100.0, 200.0, 300.0])
            assert np.allclose(ir[:, 1], [3.684, 116.334, 116.330])
            assert vcd.shape == (1, 2)

    def test_parse_freq_not_found(self, parser):
        parser.fl = "no frequencies"
        freq, ir, vcd = parser.parse_freq()
        assert len(freq) == 0
        assert ir.shape == (1, 2)
        assert vcd.shape == (1, 2)

    def test_parse_tddft_uv(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=NWCHEM_OPT + "\n" + NWCHEM_TDDFT)):
            p = NWChemParser("dummy.log", mock_logger)
            uv, ecd = p.parse_tddft()
            assert uv.shape == (2, 2)
            assert np.isclose(uv[0, 0], 5.0)
            assert np.isclose(uv[0, 1], 0.1234)
            assert ecd.shape == (1, 2)

    def test_parse_tddft_not_found(self, parser):
        parser.fl = "no excitations"
        uv, ecd = parser.parse_tddft()
        assert uv.shape == (1, 2)
        assert ecd.shape == (1, 2)
