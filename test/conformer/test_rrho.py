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

    def test_calc_S_R_grimme_fallback_and_averaging(self):
        T = 298.15
        freq = np.array([50.0, 100.0])
        # None B should fallback to GRIMME_BAV
        s_none = rrho.calc_S_R_grimme(freq, T, B=None)
        assert np.all(s_none > 0)

        # Array with zero should filter out zeros safely
        s_with_zero = rrho.calc_S_R_grimme(freq, T, B=np.array([10.0, 0.0, 5.0]))
        assert np.all(s_with_zero > 0)

        # Normal 3D rotational constants
        s_3d = rrho.calc_S_R_grimme(freq, T, B=np.array([10.0, 10.0, 10.0]))
        assert np.all(s_3d > 0)

    def test_calc_vibrational_entropy_truhlar(self):
        freq = np.array([20.0, 50.0, 300.0])
        T = 298.15
        s_truhlar = rrho.calc_vibrational_entropy_truhlar(freq, T, cut_off=100.0)
        # Should raise 20 and 50 to 100
        freq_expected = np.array([100.0, 100.0, 300.0])
        s_expected = float(np.sum(rrho.calc_S_V_grimme(freq_expected, T)) * rrho.J_TO_H)
        assert np.isclose(s_truhlar, s_expected)

    def test_rotational_entropy_symmetry_number(self):
        T = 298.15
        B = np.array([10.0, 10.0, 10.0])
        s_rot_1 = rrho.calc_rotational_entropy(B, T, symno=1)
        s_rot_2 = rrho.calc_rotational_entropy(B, T, symno=2)
        # S_rot(sigma=2) = S_rot(sigma=1) - k_B * ln(2) * J_TO_H (in Hartree/K)
        expected_diff = rrho.Boltzmann * np.log(2.0) * rrho.J_TO_H
        assert s_rot_1 - s_rot_2 == pytest.approx(expected_diff, rel=1e-5)

    def test_free_gibbs_energy_models(self):
        scf = -100.0
        T = 298.15
        freq = np.array([30.0, 150.0, 500.0])
        mw = 28.0
        B = np.array([5.0, 5.0, 5.0])
        m = 1

        G_grimme, zpve_g, h_g, S_g = rrho.free_gibbs_energy(scf, T, freq, mw, B, m, model="grimme")
        G_truhlar, zpve_t, h_t, S_t = rrho.free_gibbs_energy(scf, T, freq, mw, B, m, model="truhlar")

        assert G_grimme == pytest.approx(scf + h_g - T * S_g)
        assert G_truhlar == pytest.approx(scf + h_t - T * S_t)

        with pytest.raises(ValueError, match="Unknown quasi-RRHO model"):
            rrho.free_gibbs_energy(scf, T, freq, mw, B, m, model="unknown")