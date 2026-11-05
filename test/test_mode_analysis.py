"""Tests for the rigid-body projection of the Hessian (mode_analysis.py)."""

import numpy as np
import pytest
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations

from ensemble_analyzer.mode_analysis import project_hessian, rigid_body_basis


def _relaxed(name, fmax):
    atoms = molecule(name)
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=fmax, steps=500)
    return atoms


def _hessian(atoms, tmp_path, tag):
    vib = Vibrations(atoms, name=str(tmp_path / tag))
    vib.run()
    H = vib.get_vibrations().get_hessian_2d()
    raw = vib.get_frequencies()
    vib.clean()
    raw = np.where(np.abs(raw.imag) > 1e-8, -raw.imag, raw.real)
    return H, raw


@pytest.fixture(scope="module")
def ethanol(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ethanol")
    atoms = _relaxed("CH3CH2OH", fmax=0.05)
    H, raw = _hessian(atoms, tmp, "vib")
    return atoms, H, raw


@pytest.fixture(scope="module")
def co2(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("co2")
    atoms = _relaxed("CO2", fmax=0.001)
    H, raw = _hessian(atoms, tmp, "vib")
    return atoms, H, raw


def _rotation(seed=0):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


class TestRigidBodyBasis:
    def test_nonlinear_has_six_vectors(self, ethanol):
        atoms, _, _ = ethanol
        D = rigid_body_basis(atoms.get_masses(), atoms.get_positions())
        assert D.shape == (3 * len(atoms), 6)
        np.testing.assert_allclose(D.T @ D, np.eye(6), atol=1e-10)

    def test_linear_has_five_vectors(self, co2):
        atoms, _, _ = co2
        D = rigid_body_basis(atoms.get_masses(), atoms.get_positions())
        assert D.shape == (9, 5)

    def test_single_atom_has_three_vectors(self):
        D = rigid_body_basis(np.array([12.0]), np.zeros((1, 3)))
        assert D.shape == (3, 3)


class TestProjectHessian:
    def test_nonlinear_count(self, ethanol):
        atoms, H, _ = ethanol
        freqs, modes = project_hessian(H, atoms.get_masses(), atoms.get_positions())
        n = len(atoms)
        assert freqs.shape == (3 * n - 6,)
        assert modes.shape == (3 * n - 6, n, 3)

    def test_linear_count(self, co2):
        atoms, H, _ = co2
        freqs, modes = project_hessian(H, atoms.get_masses(), atoms.get_positions())
        assert freqs.shape == (4,)
        assert modes.shape == (4, 3, 3)

    def test_single_atom_has_no_modes(self):
        freqs, modes = project_hessian(np.zeros((3, 3)), np.array([1.0]), np.zeros((1, 3)))
        assert freqs.shape == (0,)
        assert modes.shape == (0, 1, 3)

    def test_high_frequencies_match_ase(self, ethanol):
        atoms, H, raw = ethanol
        freqs, _ = project_hessian(H, atoms.get_masses(), atoms.get_positions())
        np.testing.assert_allclose(np.sort(freqs)[-5:], np.sort(raw)[-5:], atol=0.05)

    def test_modes_orthonormal_and_internal(self, ethanol):
        atoms, H, _ = ethanol
        masses, pos = atoms.get_masses(), atoms.get_positions()
        _, modes = project_hessian(H, masses, pos)
        M = modes.reshape(len(modes), -1)
        np.testing.assert_allclose(M @ M.T, np.eye(len(modes)), atol=1e-10)
        D = rigid_body_basis(masses, pos)
        np.testing.assert_allclose(M @ D, 0.0, atol=1e-10)

    def test_invariant_under_rigid_motion(self, ethanol):
        atoms, H, _ = ethanol
        masses, pos = atoms.get_masses(), atoms.get_positions()
        ref, _ = project_hessian(H, masses, pos)

        R = _rotation(1)
        Rbig = np.kron(np.eye(len(atoms)), R)
        moved = pos @ R.T + np.array([3.0, -7.5, 12.0])
        freqs, _ = project_hessian(Rbig @ H @ Rbig.T, masses, moved)
        np.testing.assert_allclose(freqs, ref, atol=1e-6)

    def test_rigid_body_stiffness_is_removed(self, ethanol):
        """Stiffness along translations/rotations must not leak into the result."""
        atoms, H, _ = ethanol
        masses, pos = atoms.get_masses(), atoms.get_positions()
        ref, _ = project_hessian(H, masses, pos)

        D = rigid_body_basis(masses, pos)
        sqm = np.repeat(np.sqrt(masses), 3)
        # A large spurious stiffness along the rigid-body directions,
        # expressed back in the Cartesian (non mass-weighted) frame
        spurious = (sqm[:, None] * (D @ D.T) * sqm[None, :]) * 5.0
        freqs, _ = project_hessian(H + spurious, masses, pos)
        np.testing.assert_allclose(freqs, ref, atol=1e-6)

    def test_spurious_rigid_modes_present_in_raw_ase(self, ethanol):
        """ASE keeps 3N modes; the projection drops exactly the rigid ones."""
        atoms, H, raw = ethanol
        freqs, _ = project_hessian(H, atoms.get_masses(), atoms.get_positions())
        assert len(raw) - len(freqs) == 6

    def test_negative_eigenvalue_gives_negative_frequency(self):
        """A saddle along a stretch must come out as a negative frequency."""
        masses = np.array([1.0, 1.0])
        pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        k = -0.5  # eV/Å², negative curvature along the bond
        H = np.zeros((6, 6))
        e = np.array([-1.0, 0, 0, 1.0, 0, 0])
        H += k * np.outer(e, e)
        freqs, modes = project_hessian(H, masses, pos)
        assert freqs.shape == (1,)
        assert freqs[0] < 0
