"""Gaussian normal modes / displacement vectors (parser shared with FROGChem)."""
# Synthetic Gaussian freq output: water, standard orientation = rotated + translated input
# orientation; the normal modes are printed (2 decimals) in the standard orientation, as Gaussian does.
import numpy as np

_th = np.deg2rad(37.0)
ROT = np.array([[np.cos(_th), -np.sin(_th), 0.0], [np.sin(_th), np.cos(_th), 0.0], [0.0, 0.0, 1.0]]) @ \
      np.array([[1.0, 0.0, 0.0], [0.0, np.cos(1.1), -np.sin(1.1)], [0.0, np.sin(1.1), np.cos(1.1)]])
INPUT = np.array([[0.30, -0.20, 0.10], [1.20, 0.25, 0.05], [-0.35, 0.55, -0.15]])
STANDARD = (INPUT - INPUT.mean(axis=0)) @ ROT + np.array([0.0, 0.0, 0.1])
FREQ = np.array([-1638.50, 3832.69, 3945.08, 4100.00])
# true displacement vectors in the input frame: (n_modes, n_atoms, 3)
DISP_INPUT = np.array([
    [[0.05, 0.03, -0.02], [-0.40, -0.30, 0.20], [-0.39, -0.20, 0.10]],
    [[0.00, -0.04, 0.03], [0.55, 0.21, -0.11], [-0.50, 0.40, -0.30]],
    [[-0.06, 0.02, 0.01], [0.41, -0.33, 0.12], [0.52, 0.04, -0.27]],
    [[0.01, 0.07, -0.03], [-0.22, -0.48, 0.19], [0.13, -0.61, 0.36]],
])
DISP_STD = DISP_INPUT @ ROT


def _orientation(title, xyz):
    rows = '\n'.join(f'{i + 1:7d}{z:11d}{0:12d}    {x:12.6f}{y:12.6f}{w:12.6f}'
                     for i, (z, (x, y, w)) in enumerate(zip([8, 1, 1], xyz)))
    dash = ' ' + '-' * 69
    return (f'                          {title}                          \n{dash}\n'
            ' Center     Atomic      Atomic             Coordinates (Angstroms)\n'
            ' Number     Number       Type             X           Y           Z\n'
            f'{dash}\n{rows}\n{dash}\n')


def _modes(rot_str):
    out = []
    for start in range(0, len(FREQ), 3):
        idx = range(start, min(start + 3, len(FREQ)))
        out.append('                      ' + ''.join(f'{i + 1:>23d}' for i in idx))
        out.append(' Frequencies --' + ''.join(f'{FREQ[i]:>14.4f}' for i in idx))
        out.append(' Red. masses --' + ''.join(f'{1.0:>14.4f}' for i in idx))
        out.append(' IR Inten    --' + ''.join(f'{10.0 * (i + 1):>14.4f}' for i in idx))
        if rot_str:
            out.append(' Rot. str.   --' + ''.join(f'{0.5 * (i + 1):>14.4f}' for i in idx))
        out.append('  Atom  AN' + '      X      Y      Z  ' * len(idx))
        for a, z in enumerate([8, 1, 1]):
            out.append(f'{a + 1:6d}{z:4d}' + ''.join(
                ''.join(f'{v:7.2f}' for v in DISP_STD[i, a]) + '  ' for i in idx))
    return '\n'.join(out)


def make_log(rot_str=False, nosymm=False):
    head = _orientation('Input orientation:', INPUT)
    if not nosymm:
        head += _orientation('Standard orientation:', STANDARD)
    return (' Charge =  0 Multiplicity = 1\n' + head +
            ' SCF Done:  E(RB3LYP) =  -76.4089     A.U. after   10 cycles\n'
            ' Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n'
            ' activities (A**4/AMU), depolarization ratios for plane and unpolarized\n'
            ' incident light, reduced masses (AMU), force constants (mDyne/A),\n'
            ' and normal coordinates:\n' + _modes(rot_str) +
            '\n\n\n -------------------\n - Thermochemistry -\n -------------------\n'
            ' Normal termination of Gaussian 16\n')


# ---------------------------------------------------------------------------
import pytest
from unittest.mock import MagicMock, mock_open, patch
from ensemble_analyzer.parsers.gaussian import GaussianParser


def _parser(text):
    with patch("builtins.open", mock_open(read_data=text)):
        return GaussianParser("dummy.log", MagicMock())


def test_displacements_rotated_to_input_orientation():
    p = _parser(make_log())
    freq, _, _ = p.parse_freq()
    modes = p.parse_normal_modes(n_atoms=3)
    assert np.allclose(p.parse_geom(), INPUT)
    assert modes.shape == (len(freq), 3, 3) == DISP_INPUT.shape
    # printed with 2 decimals in the standard orientation -> rounding error only
    assert np.allclose(modes, DISP_INPUT, atol=0.01)
    assert not np.allclose(modes, DISP_STD, atol=0.05)


def test_nosymm_modes_unchanged():
    modes = _parser(make_log(nosymm=True)).parse_normal_modes(n_atoms=3)
    assert np.allclose(modes, np.round(DISP_STD, 2))


def test_freq_without_rot_str():
    freq, ir, vcd = _parser(make_log(rot_str=False)).parse_freq()
    assert np.allclose(freq, FREQ)
    assert np.allclose(ir[:, 1], [10.0, 20.0, 30.0, 40.0])
    assert vcd.shape == (1, 2) and not vcd.any()


def test_freq_with_rot_str():
    freq, ir, vcd = _parser(make_log(rot_str=True)).parse_freq()
    assert vcd.shape == (4, 2)
    assert np.allclose(vcd[:, 1], [0.5, 1.0, 1.5, 2.0])
    # normal modes are still read when the VCD rows are present
    assert np.allclose(_parser(make_log(rot_str=True)).parse_normal_modes(3), DISP_INPUT, atol=0.01)
