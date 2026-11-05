"""BaseMlCalc._run_vibrations: rigid-body modes must be projected out."""

from types import SimpleNamespace

import numpy as np
import pytest
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.emt import EMT
from ase.optimize import BFGS

from ensemble_analyzer.calculators._ml_base import BaseMlCalc
from ensemble_analyzer.mode_analysis import NormalModeAnalyzer


class _EMTWithDipole(Calculator):
    """EMT plus a point-charge dipole, so that ase.vibrations.Infrared works."""

    implemented_properties = ["energy", "forces", "dipole"]

    def __init__(self, charges, **kwargs):
        super().__init__(**kwargs)
        self._emt = EMT()
        self._charges = np.asarray(charges, dtype=float)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        a = self.atoms.copy()
        a.calc = self._emt
        self.results["energy"] = a.get_potential_energy()
        self.results["forces"] = a.get_forces()
        self.results["dipole"] = self._charges @ self.atoms.get_positions()


def _make_calc(tmp_path):
    protocol = SimpleNamespace(number=3, constraints=None)
    conf = SimpleNamespace(folder="conf_1")
    (tmp_path / "conf_1").mkdir()
    return BaseMlCalc(protocol=protocol, cpu=1, conf=conf)


def _relaxed_ethanol(calc):
    atoms = molecule("CH3CH2OH")
    atoms.calc = calc
    BFGS(atoms, logfile=None).run(fmax=0.05, steps=500)
    return atoms


def test_vibrations_without_dipole_are_projected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = _relaxed_ethanol(EMT())
    freqs, intensities, modes = _make_calc(tmp_path)._run_vibrations(atoms)

    n_modes = 3 * len(atoms) - 6
    assert freqs.shape == (n_modes,)
    assert intensities.shape == (n_modes,)
    assert np.all(intensities == 0)
    assert modes.shape == (n_modes, len(atoms), 3)
    # NormalModeAnalyzer drops zero-norm modes: none must be dropped,
    # otherwise mode indices would no longer match Freq
    analyzer = NormalModeAnalyzer(modes, atoms.get_positions(), tuple(atoms.get_chemical_symbols()))
    assert analyzer.n_modes == n_modes


def test_vibrations_with_dipole_map_intensities(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    symbols = molecule("CH3CH2OH").get_chemical_symbols()
    charges = [{"C": -0.2, "O": -0.6, "H": 0.0}[s] for s in symbols]
    # make it neutral by spreading the compensating charge on the hydrogens
    n_h = symbols.count("H")
    charges = [c + (0.8 / n_h if s == "H" else 0.0) for c, s in zip(charges, symbols)]
    atoms = _relaxed_ethanol(_EMTWithDipole(charges))

    freqs, intensities, modes = _make_calc(tmp_path)._run_vibrations(atoms)

    n_modes = 3 * len(atoms) - 6
    assert freqs.shape == intensities.shape == (n_modes,)
    assert np.all(intensities >= 0)
    assert intensities.sum() > 0
    cache = tmp_path / "conf_1" / "protocol_3"
    assert not cache.exists() or not any(cache.iterdir())  # cache cleaned


@pytest.mark.parametrize("seed", [0, 1])
def test_frequencies_invariant_to_rigid_motion(tmp_path, monkeypatch, seed):
    monkeypatch.chdir(tmp_path)
    atoms = _relaxed_ethanol(EMT())
    ref, _, _ = _make_calc(tmp_path)._run_vibrations(atoms)

    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    moved = atoms.copy()
    moved.positions = atoms.positions @ q.T + rng.normal(size=3) * 5
    moved.calc = EMT()
    (tmp_path / "conf_1").rename(tmp_path / f"old_{seed}")
    freqs, _, _ = _make_calc(tmp_path)._run_vibrations(moved)

    # The finite-difference Hessian itself depends slightly on orientation
    # (Cartesian displacement directions), so soft modes move by a few cm⁻¹;
    # exact invariance of the projection is tested in test_mode_analysis.py.
    assert freqs.shape == ref.shape
    stiff = np.abs(ref) > 200
    np.testing.assert_allclose(freqs[stiff], ref[stiff], atol=1.0)
