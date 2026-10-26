"""Tests for the optional PRISM Pruner backend and the eigenvalue-ordering fix."""

import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from ensemble_analyzer._managers import prism_backend
from ensemble_analyzer._managers.prism_backend import kabsch_rmsd, run_prism_pruning
from ensemble_analyzer._managers.pruning_manager import PruningManager
from ensemble_analyzer.constants import EH_TO_KCAL

import importlib.util

HAS_PRISM = importlib.util.find_spec("prism_pruner") is not None
needs_prism = pytest.mark.skipif(not HAS_PRISM, reason="prism_pruner not installed")


# ---------------------------------------------------------------------------
# Small chiral test system (ethanol with random distortion -> C1 symmetry)
# ---------------------------------------------------------------------------

ATOMS = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]


@pytest.fixture
def structures():
    rng = np.random.default_rng(0)
    base = np.array([
        [ 1.168, -0.400, 0.000], [ 0.000,  0.559, 0.000], [-1.190, -0.229, 0.000],
        [-1.946,  0.371, 0.000], [ 0.042,  1.201, 0.880], [ 0.042,  1.201, -0.880],
        [ 2.115,  0.150, 0.000], [ 1.128, -1.040, 0.885], [ 1.128, -1.040, -0.885],
    ])
    X = base + rng.normal(scale=0.08, size=base.shape)
    B = X.copy()
    B[0] += [1.5, 0.0, 0.0]
    B[-1] += [1.2, 0.6, 0.0]
    X_dup = X + rng.normal(scale=0.003, size=X.shape)
    mirror = X.copy()
    mirror[:, 0] *= -1
    return SimpleNamespace(X=X, B=B, X_dup=X_dup, mirror=mirror)


class FakeConf:
    """Minimal stand-in for Conformer as used by PruningManager."""

    def __init__(self, number, geom, energy_kcal):
        self.number = number
        self.atoms = tuple(ATOMS)
        self.last_geometry = np.asarray(geom, float)
        self.active = True
        self.diactivated_by = None
        self.rotatory = 1.0
        self.moment = 1.0
        e = energy_kcal / EH_TO_KCAL
        self.energies = SimpleNamespace(get_energy=lambda e=e: e)

    def distance_matrix(self, include_H, geom=None):
        g = self.last_geometry
        return np.linalg.norm(g[:, None] - g[None], axis=-1)


def prism_protocol(**kw):
    base = dict(prune_method="prism", thrG=1.0, number=1, no_prune=False,
                prism_moi=False, prism_moi_dev=1e-2, prism_max_rmsd=0.25,
                prism_max_dev=None, prism_timeout=60)
    base.update(kw)
    return SimpleNamespace(**base)


# ---------------------------------------------------------------------------
# Backend (requires prism_pruner)
# ---------------------------------------------------------------------------

@needs_prism
@pytest.mark.parametrize("method,use_moi", [("prism", True), ("prism", False), ("prism_rotcorr", True)])
def test_keeps_lowest_energy_duplicate(structures, method, use_moi):
    s = structures
    coords = np.array([s.X_dup, s.X, s.B])          # X_dup is lowest in energy
    keep, records = run_prism_pruning(coords, ATOMS, np.array([0.0, 0.2, 0.5]), max_dE=1.0,
                                      method=method, use_moi=use_moi, heavy_atoms_only=False)
    assert keep.tolist() == [True, False, True]
    assert [(r.removed, r.reference) for r in records] == [(1, 0)]


@needs_prism
def test_unique_structure_never_lost(structures):
    """Regression for the prism_pruner 0.1.2 internal re-sorting issue."""
    s = structures
    rng = np.random.default_rng(1)
    for _ in range(50):
        e = np.sort(rng.uniform(0, 0.9, 3))
        perm = rng.permutation(3)
        coords = np.array([s.X, s.B, s.X_dup])[perm]
        keep, _ = run_prism_pruning(coords, ATOMS, e, max_dE=1.0, heavy_atoms_only=False)
        kept = coords[keep]
        assert any(np.allclose(k, s.B) for k in kept)
        assert keep.sum() == 2
        dup_pos = [i for i in range(3) if not np.allclose(coords[i], s.B)]
        assert keep[min(dup_pos)]  # lowest-energy duplicate retained


@needs_prism
def test_energy_window_blocks_comparison(structures):
    s = structures
    keep, records = run_prism_pruning(np.array([s.X, s.X_dup]), ATOMS, np.array([0.0, 2.0]),
                                      max_dE=1.0, heavy_atoms_only=False)
    assert keep.all() and records == []


@needs_prism
def test_mirror_image_handling(structures):
    s = structures
    coords = np.array([s.X, s.mirror])
    e = np.array([0.0, 0.0])
    keep_moi, _ = run_prism_pruning(coords, ATOMS, e, max_dE=1.0, use_moi=True, heavy_atoms_only=False)
    keep_rmsd, _ = run_prism_pruning(coords, ATOMS, e, max_dE=1.0, use_moi=False, heavy_atoms_only=False)
    assert keep_moi.sum() == 1          # MOI is reflection-invariant: mirror merged
    assert keep_rmsd.all()              # Kabsch RMSD (proper rotations) keeps both


@needs_prism
def test_rejects_unsorted_input(structures):
    s = structures
    with pytest.raises(ValueError, match="ascending"):
        run_prism_pruning(np.array([s.X, s.B]), ATOMS, np.array([1.0, 0.0]), max_dE=1.0)


@needs_prism
def test_manager_prism_integration(mock_logger, structures):
    s = structures
    # Deliberately unsorted list: the manager must sort by energy itself.
    confs = [FakeConf(1, s.X, 0.2), FakeConf(2, s.B, 0.5), FakeConf(3, s.X_dup, 0.0)]
    mgr = PruningManager(logger=mock_logger, include_H=True)
    mgr._remove_duplicates(confs, prism_protocol())
    assert [c.active for c in confs] == [False, True, True]
    assert confs[0].diactivated_by == 3
    rec = mgr._deactivation_records
    assert len(rec) == 1 and rec[0].check_id == 1 and rec[0].reference_id == 3
    assert rec[0].criterion.startswith("prism:")
    assert "Criterion" in rec[0].to_dict()


# ---------------------------------------------------------------------------
# Always-run tests (no prism needed)
# ---------------------------------------------------------------------------

def test_dispatch_defaults_to_rotconst(mock_logger):
    mgr = PruningManager(logger=mock_logger)
    with patch.object(mgr, "_remove_duplicates_rotconst") as rc, \
         patch.object(mgr, "_remove_duplicates_prism") as pr:
        mgr._remove_duplicates([], SimpleNamespace(thrG=1, thrB=1, number=1))   # no attribute
        mgr._remove_duplicates([], MagicMock())                                 # mock attribute
        mgr._remove_duplicates([], SimpleNamespace(prune_method="rotconst"))
    assert rc.call_count == 3 and pr.call_count == 0


def test_dispatch_to_prism(mock_logger):
    mgr = PruningManager(logger=mock_logger)
    with patch.object(mgr, "_remove_duplicates_prism") as pr:
        mgr._remove_duplicates([], SimpleNamespace(prune_method="PRISM_RotCorr"))
    pr.assert_called_once()
    assert pr.call_args.args[2] == "prism_rotcorr"


def test_import_error_message(monkeypatch):
    monkeypatch.setitem(sys.modules, "prism_pruner", None)
    with pytest.raises(ImportError, match=r"pip install prism_pruner"):
        prism_backend.import_prism()


def test_invalid_method_rejected():
    with pytest.raises(ValueError):
        run_prism_pruning(np.zeros((2, 3, 3)), ["C"] * 3, np.zeros(2), max_dE=1.0, method="foo")


def test_protocol_validates_prune_method():
    from ensemble_analyzer.protocol.protocol import Protocol
    p = Protocol(number=0, functional="GFN2-xTB", calculator="tblite", prune_method=" Prism ")
    assert p.prune_method == "prism"
    assert p.prism_moi is False          # MOI stage is opt-in
    with pytest.raises(ValueError, match="prune_method"):
        Protocol(number=0, functional="GFN2-xTB", calculator="tblite", prune_method="bogus")
    with pytest.raises(ValueError):
        Protocol(number=0, functional="GFN2-xTB", calculator="tblite", prism_max_rmsd=0)


def test_kabsch_rmsd_rotation_invariant_and_mirror_sensitive(structures):
    s = structures
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    assert kabsch_rmsd(s.X, s.X @ R.T + 3.0) < 1e-8
    assert kabsch_rmsd(s.X, s.mirror) > 0.1


# ---------------------------------------------------------------------------
# Eigenvalue-ordering fix (pruning log metric + clustering features)
# ---------------------------------------------------------------------------

def test_lambda_rmsd_small_for_near_duplicates():
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(100):
        x = rng.normal(size=(20, 3)) * 2
        a, b = FakeConf(1, x, 0.0), FakeConf(2, x + rng.normal(scale=0.01, size=x.shape), 0.0)
        a.atoms = b.atoms = tuple(["C"] * 20)
        worst = max(worst, PruningManager._calculate_rmsd(a, b, include_H=True))
    assert worst < 0.1


def test_cluster_features_sorted_and_permutation_invariant(mock_logger):
    from ensemble_analyzer._clustering.cluster_manager import ClusteringManager
    rng = np.random.default_rng(0)
    x = rng.normal(size=(12, 3))
    perm = rng.permutation(12)
    geoms = np.array([x, x[perm]])
    atoms = np.array(["C"] * 12)
    mgr = ClusteringManager.__new__(ClusteringManager)
    feats = mgr._calculate_distance_matrix_eigenvalues(geoms, atoms, include_H=True)
    assert np.all(np.diff(feats, axis=1) >= -1e-12)
    np.testing.assert_allclose(feats[0], feats[1], atol=1e-10)
