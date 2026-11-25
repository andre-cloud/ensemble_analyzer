"""Point group and NMR equivalence tests, ported from FROGChem (test_symmetry_enhanced.py, test_nmr_sim.py)."""

import numpy as np
import pytest

from ensemble_analyzer._nmr.equivalence import PointGroupDetector, detect_point_group, get_nmr_equivalencies


class TestPointGroupDetector:
    def test_water_c2v(self):
        # H2O in yz plane: O at origin, H at (+-0.757, 0.586, 0)
        symbols = ['O', 'H', 'H']
        coords = np.array([
            [0.0, 0.0, 0.117],
            [0.0, 0.757, -0.469],
            [0.0, -0.757, -0.469]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['is_valid'] is True
        assert res['sch_symbol'] == 'C2v'
        assert 'H' in res['equivalencies']
        # Both H atoms should be in the same equivalent set
        assert len(res['equivalencies']['H']) == 1
        assert sorted(res['equivalencies']['H'][0]) == [1, 2]

    def test_ammonia_c3v(self):
        # NH3: C3 axis along z
        symbols = ['N', 'H', 'H', 'H']
        coords = np.array([
            [0.0, 0.0, 0.116],
            [0.0, 0.938, -0.271],
            [0.812, -0.469, -0.271],
            [-0.812, -0.469, -0.271]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'C3v'
        assert len(res['equivalencies']['H']) == 1
        assert sorted(res['equivalencies']['H'][0]) == [1, 2, 3]

    def test_planar_d3h(self):
        # BF3: Boron at origin, 3 F atoms in xy plane 120 deg apart
        symbols = ['B', 'F', 'F', 'F']
        r = 1.31
        coords = np.array([
            [0.0, 0.0, 0.0],
            [r, 0.0, 0.0],
            [-0.5 * r, np.sqrt(3)/2 * r, 0.0],
            [-0.5 * r, -np.sqrt(3)/2 * r, 0.0]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'D3h'

    def test_linear_diatomic(self):
        # CO: C_inf_v
        symbols = ['C', 'O']
        coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.13]])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'C_inf_v'

    def test_linear_centrosymmetric(self):
        # CO2: D_inf_h
        symbols = ['O', 'C', 'O']
        coords = np.array([[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'D_inf_h'

    def test_methane_td(self):
        # CH4: Td
        symbols = ['C', 'H', 'H', 'H', 'H']
        a = 1.09 / np.sqrt(3)
        coords = np.array([
            [0.0, 0.0, 0.0],
            [a, a, a],
            [a, -a, -a],
            [-a, a, -a],
            [-a, -a, a]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'Td'

    def test_asymmetric_c1(self):
        # CHFClBr
        symbols = ['C', 'H', 'F', 'Cl', 'Br']
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.09],
            [1.35, 0.0, -0.4],
            [-0.8, 1.5, -0.4],
            [-0.7, -1.6, -0.5]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'C1'

    def test_planar_cs(self):
        # Hypochlorous acid HOCl (planar, no other symmetry => Cs)
        symbols = ['O', 'H', 'Cl']
        coords = np.array([
            [0.0, 0.0, 0.0],
            [0.96, 0.0, 0.0],
            [-0.5, 1.6, 0.0]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'Cs'

    def test_c2h(self):
        # trans-1,2-dichloroethene-like symmetry: C2h
        symbols = ['C', 'C', 'Cl', 'Cl']
        coords = np.array([
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [1.0, 1.5, 0.5],
            [-1.0, -1.5, -0.5]
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'C2h'

    def test_pure_ci(self):
        # 3 pairs of distinct atoms related ONLY by inversion
        symbols = ['F', 'F', 'Cl', 'Cl', 'Br', 'Br']
        coords = np.array([
            [1.2, 2.3, 0.7],
            [-1.2, -2.3, -0.7],
            [2.1, -1.5, 3.2],
            [-2.1, 1.5, -3.2],
            [-3.4, 0.8, 1.9],
            [3.4, -0.8, -1.9],
        ])
        res = detect_point_group(symbols, coords, tolerance=0.1, exclude_h=False)
        assert res['sch_symbol'] == 'Ci'


def test_nmr_equivalencies_methane():
    # Tetrahedral methane
    atoms = ['C', 'H', 'H', 'H', 'H']
    coords = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.09],
        [1.027, 0.0, -0.363],
        [-0.513, 0.890, -0.363],
        [-0.513, -0.890, -0.363],
    ]
    equiv = get_nmr_equivalencies(atoms, coords)
    assert 'H' in equiv
    # All 4 hydrogens must be grouped together into 1 equivalent group of 4
    assert len(equiv['H']) == 1
    assert set(equiv['H'][0]) == {1, 2, 3, 4}


def test_nmr_equivalencies_methyl_rotor():
    # Molecule with a methyl group: C1-C2 where C2 is bonded to 3 H
    atoms = ['C', 'C', 'H', 'H', 'H', 'H']
    coords = [
        [0.0, 0.0, 0.0],
        [1.54, 0.0, 0.0],
        [-0.5, 0.9, 0.0],      # H on C1
        [1.9, 0.9, 0.0],       # H on C2
        [1.9, -0.45, 0.78],    # H on C2
        [1.9, -0.45, -0.78],   # H on C2
    ]
    equiv = get_nmr_equivalencies(atoms, coords)
    assert 'H' in equiv
    # The 3 methyl protons (indices 3, 4, 5) must be grouped together
    methyl_found = False
    for group in equiv['H']:
        if set(group) == {3, 4, 5}:
            methyl_found = True
    assert methyl_found, f"Methyl rotor grouping failed, got: {equiv['H']}"


def test_equivalence_rejects_heavy_atom_symmetry_broken_by_h():
    atoms = ['C', 'C', 'H', 'H', 'H', 'H']
    coords = [[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [-0.5, 0.9, 0.0],
              [1.9, 0.9, 0.0], [1.9, -0.45, 0.78], [1.9, -0.45, -0.78]]
    equiv = get_nmr_equivalencies(atoms, coords)
    groups = [set(g) for g in equiv['H']]
    assert {3, 4, 5} in groups and {2} in groups
    assert [set(g) for g in equiv['C']] == [{0}, {1}]
