"""
Tests for the Gaussian Parser module.
Verifies regex extraction of geometry, frequencies, energy, and TD-DFT data.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, mock_open, patch
from ensemble_analyzer._parsers._gaussian import GaussianParser

class TestGaussianParser:
    
    @pytest.fixture
    def parser(self, mock_logger):
        """Fixture providing an instance of GaussianParser with mocked file IO."""
        with patch("builtins.open", mock_open(read_data="")):
             p = GaussianParser("dummy.log", mock_logger)
             p.fl = ""
             return p

    def test_parse_geom_regex(self, parser):
        """Test extraction of the final geometry."""
        parser.fl = """
 Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        1.000000    2.000000    3.000000
 ---------------------------------------------------------------------
 --"""
        geom = parser.parse_geom()
        assert np.allclose(geom, [[1.0, 2.0, 3.0]])

    def test_parse_freq_parsing(self, parser):
        """Test extraction of frequency, IR intensity, and VCD (Rot. str.) data."""
        parser.fl = """
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 Frequencies --   100.0000   200.0000
 IR Inten    --     5.0000    10.0000
 Rot. str.   --     0.1000    -0.2000
 
 
 
 """
        freq, ir, vcd = parser.parse_freq()
        assert np.allclose(freq, [100.0, 200.0])
        assert ir.shape == (2, 2) # [freq, inten]
        assert vcd.shape == (2, 2) # [freq, rot]
        assert ir[0, 1] == 5.0

    def test_parse_tddft_missing(self, parser):
        """Test behavior when TD-DFT data is missing."""
        parser.fl = "No excited states here"
        uv, ecd = parser.parse_tddft()
        assert np.array_equal(uv, np.zeros((1, 2)))

    def test_parse_energy(self, parser):
        """Test extraction of SCF Energy."""
        parser.fl = """
        Step number   10       Rms gradient =  0.12345678
        SCF Done:  E(RB3LYP) =  -1234.56789012     A.U. after   11 cycles
        """
        assert parser.parse_energy() == -1234.56789012

    def test_parse_B_m(self, parser, mock_logger):
        """Test extraction and conversion of Rotational Constants and Dipole Moment."""
        parser.fl = """
        Rotational constants (GHZ):      10.00000      5.00000      2.00000
        Standard basis: 6-31G(d)
        Dipole moment (field-independent basis, Debye):
        X=     1.5000    Y=     -0.5000    Z=     0.0000  Tot=     1.5811
        Quadrupole moment (field-independent basis, Debye-Ang):
        """
        
        # Patch constants to ensure deterministic conversion check
        with patch("ensemble_analyzer._parsers._gaussian.CONVERT_B", {"GHz": 29.9792458}):
            B, M = parser.parse_B_m()
        
        # Verify Dipole
        assert np.allclose(M, [1.5, -0.5, 0.0])
        
        # Verify Rotational Constants (10.0 / factor)
        expected_B0 = 10.0 / 29.9792458
        assert np.isclose(B[0], expected_B0)

    def test_parse_tddft_found(self, parser):
        """Test parsing of UV and ECD spectra from TD-DFT output."""
        parser.fl = """
        Excited states from <AA,BB,CC,DD> symmetry
        
        Excitation energies and oscillator strengths:
        Excited State   1:      Singlet-A      3.5000 eV  354.24 nm  f=0.1234 <S**2>=0.000
        Excited State   2:      Singlet-A      4.0000 eV  309.96 nm  f=0.5678 <S**2>=0.000
        ...
        ***
        
        1/2[<0|del|b> * <b|rxdel|0> + <0|del|b> * <b|delr+rdel|0>]
         State         X           Y           Z        R(velocity)
           1         0.1111     -0.2222      0.3333     -12.3456
           2         0.4444      0.5555      0.6666      78.9012
         1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]
        """
        
        uv, ecd = parser.parse_tddft()
        
        assert uv.shape == (2, 2)
        assert uv[0, 0] == 3.5    # Energy
        assert uv[0, 1] == 0.1234 # Oscillator strength
        
        assert ecd.shape == (2, 2)
        assert ecd[0, 1] == -12.3456 # Rotatory strength