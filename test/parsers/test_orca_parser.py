import pytest
import numpy as np
from unittest.mock import MagicMock, mock_open, patch
from ensemble_analyzer._parsers._orca import OrcaParser


ORCA5_OPT = """
Program Version 5

ORCA TERMINATED NORMALLY
"""

ORCA6_OPT = """
Program Version 6

ORCA TERMINATED NORMALLY
"""

ORCA6_OPT_CONVERGED = """
Program Version 6

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  C    0.000000    0.000000    0.000000
  O    0.000000    0.000000    1.200000

THE OPTIMIZATION HAS CONVERGED

ORCA TERMINATED NORMALLY
"""

ORCA6_ENERGY = """
Program Version 6

FINAL SINGLE POINT ENERGY    -113.45678900

ORCA TERMINATED NORMALLY
"""

ORCA6_BCONST = """
Program Version 6

Rotational constants in cm-1:       1.50000       0.20000       0.10000
Total Dipole Moment   :     0.50000     0.30000     0.10000

ORCA TERMINATED NORMALLY
"""

ORCA6_FREQ = """
Program Version 6

VIBRATIONAL FREQUENCIES
------------
   0:      500.00
   1:     1000.00
   2:       -50.00



Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
   0    500.00    1.00    10.00    0.0100    0.1000    0.0000    0.0000
   1   1000.00    2.00    20.00    0.0200    0.2000    0.0000    0.0000



Mode   Freq    VCD-Intensity    
       (1/cm) (1E-44*esu^2*cm^2) 
---------------------------------
   0    500.00    0.50
   1   1000.00   -0.30



ORCA TERMINATED NORMALLY
"""

ORCA5_FREQ = """
Program Version 5

VIBRATIONAL FREQUENCIES
------------
   0:      500.00
   1:     1000.00

Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
   0    500.00    1.00    10.00    0.0100    0.1000    0.0000    0.0000
   1   1000.00    2.00    20.00    0.0200    0.2000    0.0000    0.0000

ORCA TERMINATED NORMALLY
"""

ORCA6_TDDFT = """
Program Version 6

SPECTRA
ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  
----------------------------------------------------------------------------------------------------
          1          3.5000  28233.00   354.24  0.1234    0.5000    0.1000   -0.2000    0.3000
          2          4.0000  32266.00   309.96  0.5678    0.8000    0.4000    0.5000    0.6000

CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ   
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)  
------------------------------------------------------------------------------------------
          1          3.5000  28233.00   354.24  12.3456    0.1000    0.2000    0.3000
          2          4.0000  32266.00   309.96 -78.9012    0.4000    0.5000    0.6000

***
ORCA TERMINATED NORMALLY
"""

ORCA5_TDDFT = """
Program Version 5

SPECTRA
ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy    Wavelength  fosc         T2        TX        TY        TZ  
        (cm-1)      (nm)                 (au**2)    (au)      (au)      (au) 
-----------------------------------------------------------------------------
   1  28233.00     354.24    0.1234    0.5000    0.1000   -0.2000    0.3000
   2  32266.00     309.96    0.5678    0.8000    0.4000    0.5000    0.6000

CD SPECTRUM
-------------------------------------------------------------------
State  Energy     Wavelength     R         MX        MY        MZ   
       (cm-1)       (nm)     (1e40*cgs)   (au)      (au)      (au)  
-------------------------------------------------------------------
   1  28233.00     354.24     12.3456    0.1000    0.2000    0.3000
   2  32266.00     309.96    -78.9012    0.4000    0.5000    0.6000

***
ORCA TERMINATED NORMALLY
"""

ORCA_CRASHED = """
Program Version 6

Some partial output but no termination
"""


class TestOrcaParser:

    @pytest.fixture
    def parser_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            p.fl = ""
            return p

    def test_version_detection_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.version == "6"

    def test_version_detection_v5(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.version == "5"

    def test_normal_termination(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.normal_termination() is True

    def test_crashed_calculation(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA_CRASHED)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.correct_exiting is False
            assert p.normal_termination() is False
            p.log.warning.assert_called()

    def test_parse_geom(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT_CONVERGED)):
            p = OrcaParser("dummy.out", mock_logger)
            geom = p.parse_geom()
            assert geom.shape == (2, 3)
            assert np.allclose(geom[0], [0.0, 0.0, 0.0])
            assert np.allclose(geom[1], [0.0, 0.0, 1.2])

    def test_opt_done(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT_CONVERGED)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.opt_done() is True

    def test_opt_not_done(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.opt_done() is False

    def test_parse_energy_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_ENERGY)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.parse_energy() == -113.45678900

    def test_parse_energy_not_found(self, parser_v6):
        parser_v6.fl = "no energy here"
        assert parser_v6.parse_energy() == 0.0

    def test_parse_B_m(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_BCONST)):
            p = OrcaParser("dummy.out", mock_logger)
            B, M = p.parse_B_m()
            assert np.allclose(B, [1.5, 0.2, 0.1])
            assert np.allclose(M, [0.5, 0.3, 0.1])

    def test_parse_B_m_not_found(self, parser_v6):
        parser_v6.fl = "no B or M here"
        B, M = parser_v6.parse_B_m()
        assert np.allclose(B, [1, 0, 0])
        assert np.allclose(M, [1, 0, 0])

    def test_parse_freq_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_FREQ)):
            p = OrcaParser("dummy.out", mock_logger)
            freq, ir, vcd = p.parse_freq()
            assert len(freq) == 3
            assert np.allclose(freq, [500.0, 1000.0, -50.0])
            assert ir.shape[0] == 2
            assert vcd.shape[0] == 2

    def test_parse_freq_no_vcd_in_v5(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_FREQ)):
            p = OrcaParser("dummy.out", mock_logger)
            freq, ir, vcd = p.parse_freq()
            assert vcd.shape == (1, 2)

    def test_parse_freq_not_found(self, parser_v6):
        parser_v6.fl = "no frequencies"
        freq, ir, vcd = parser_v6.parse_freq()
        assert len(freq) == 0
        assert ir.shape == (1, 2)
        assert vcd.shape == (1, 2)

    def test_parse_tddft_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_TDDFT)):
            p = OrcaParser("dummy.out", mock_logger)
            uv, ecd = p.parse_tddft()
            assert uv.shape == (2, 2)
            assert np.isclose(uv[0, 0], 354.24)
            assert np.isclose(uv[0, 1], 0.1000)
            assert ecd.shape == (2, 2)

    def test_parse_tddft_v5(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_TDDFT)):
            p = OrcaParser("dummy.out", mock_logger)
            uv, ecd = p.parse_tddft()
            assert uv.shape == (2, 2)
            assert ecd.shape == (2, 2)

    def test_parse_tddft_not_found(self, parser_v6):
        parser_v6.fl = "no spectra"
        uv, ecd = parser_v6.parse_tddft()
        assert uv.shape == (1, 2)
        assert ecd.shape == (1, 2)
