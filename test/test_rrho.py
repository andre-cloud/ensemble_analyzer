import pytest
import numpy as np
from unittest.mock import patch
from ensemble_analyzer import rrho

class TestRRHO:
    
    @pytest.fixture
    def mock_constants(self):
        with patch.multiple("ensemble_analyzer.rrho", 
                            h=1.0, c=1.0, Boltzmann=1.0, J_TO_H=1.0, N_A=1.0):
            yield

    def test_calc_damp(self):
        res = rrho.calc_damp(np.array([100.0]), cut_off=100.0, alpha=1)
        assert res[0] == 0.5
        res = rrho.calc_damp(np.array([1000.0]), cut_off=1.0, alpha=1)
        assert res[0] > 0.99

    def test_calc_zpe(self, mock_constants):
        freqs = np.array([10.0, 20.0])
        zpe = rrho.calc_zpe(freqs)
        assert zpe == 15.0

    def test_calc_translational_energy(self, mock_constants):
        T = 2.0
        assert rrho.calc_translational_energy(T) == 3.0

    def test_calc_rotational_energy(self, mock_constants):
        T = 2.0
        assert rrho.calc_rotational_energy(T, linear=False) == 3.0
        assert rrho.calc_rotational_energy(T, linear=True) == 2.0

    def test_calc_qRRHO_energy(self, mock_constants):
        freq = np.array([1.0]) 
        T = 1.0
        val = rrho.calc_qRRHO_energy(freq, T)
        expected = 1.0 * np.exp(-1) / (1 - np.exp(-1))
        assert np.isclose(val[0], expected)

    def test_free_gibbs_energy_integration(self):
        scf = -100.0
        T = 298.15
        freq = np.array([100.0, 200.0, 300.0])
        mw = 18.0 
        B = np.array([10.0, 10.0, 10.0])
        m = 1
        
        G, zpve, H_corr, S = rrho.free_gibbs_energy(scf, T, freq, mw, B, m)
        
        # H_corr is the thermal correction, so H_total = SCF + H_corr
        H_total = scf + H_corr
        
        # Verify Gibbs relation: G = H_total - T*S
        assert G == pytest.approx(H_total - T * S)