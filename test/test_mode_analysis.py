import numpy as np
from ensemble_analyzer.mode_analysis import NormalModeAnalyzer


class TestDeriveInternalsFromFragments:
    def test_single_bond(self):
        frags = [[0, 1]]
        result = NormalModeAnalyzer.derive_internals_from_fragments(frags)
        assert result == [[0, 1]]

    def test_triatomic_fragment(self):
        frags = [[0, 1, 2]]
        result = NormalModeAnalyzer.derive_internals_from_fragments(frags)
        assert [0, 1] in result
        assert [1, 2] in result
        assert [0, 1, 2] in result
        assert len(result) == 3

    def test_four_atom_fragment(self):
        frags = [[0, 1, 2, 3]]
        result = NormalModeAnalyzer.derive_internals_from_fragments(frags)
        assert [0, 1] in result
        assert [1, 2] in result
        assert [2, 3] in result
        assert [0, 1, 2] in result
        assert [1, 2, 3] in result
        assert [0, 1, 2, 3] in result
        assert len(result) == 6

    def test_multiple_fragments(self):
        frags = [[0, 1, 2], [3, 4]]
        result = NormalModeAnalyzer.derive_internals_from_fragments(frags)
        assert [0, 1] in result
        assert [1, 2] in result
        assert [0, 1, 2] in result
        assert [3, 4] in result
        assert len(result) == 4

    def test_empty_fragments(self):
        result = NormalModeAnalyzer.derive_internals_from_fragments([])
        assert result == []

    def test_single_atom_fragment(self):
        result = NormalModeAnalyzer.derive_internals_from_fragments([[5]])
        assert result == []

    def test_only_bonds_for_pair(self):
        frags = [[0, 1], [2, 3]]
        result = NormalModeAnalyzer.derive_internals_from_fragments(frags)
        assert result == [[0, 1], [2, 3]]


class TestProjectOnInternals:
    def test_bond_stretch(self):
        """Bond stretching along the mode should give positive projection."""
        atoms = ("C", "O")
        geom = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
        modes = np.zeros((1, 2, 3))
        modes[0, 0] = [-0.5, 0.0, 0.0]
        modes[0, 1] = [0.5, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        result = ana.project_on_internals(0, [[0, 1]])  # bond C-O
        assert len(result) == 1
        label, comp, pct = result[0]
        assert label == "B(C0,O1)"
        assert comp > 0
        assert abs(comp - 1.0) < 0.01
        assert abs(pct - 100.0) < 0.01

    def test_bond_compression(self):
        """Bond compression along the mode should give negative projection."""
        atoms = ("C", "O")
        geom = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
        modes = np.zeros((1, 2, 3))
        modes[0, 0] = [0.5, 0.0, 0.0]
        modes[0, 1] = [-0.5, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        result = ana.project_on_internals(0, [[0, 1]])
        label, comp, pct = result[0]
        assert comp < 0
        assert abs(comp - (-1.0)) < 0.01
        assert abs(pct - 100.0) < 0.01

    def test_angle_bend(self):
        """Bending mode should project onto angle internal."""
        atoms = ("O", "C", "O")
        geom = np.array([
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.5, 0.866, 0.0],
        ])
        modes = np.zeros((1, 3, 3))
        modes[0, 0] = [0.0, 0.0, 0.0]
        modes[0, 1] = [0.0, 0.0, 0.0]
        modes[0, 2] = [0.0, 0.01, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        result = ana.project_on_internals(0, [[0, 1, 2]])
        assert len(result) == 1
        label, comp, pct = result[0]
        assert label == "A(O0,C1,O2)"
        assert abs(pct - 100.0) < 0.01

    def test_dihedral_torsion(self):
        """Torsional mode should project onto dihedral internal."""
        atoms = ("C", "C", "C", "C")
        geom = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.5, 0.5, 0.0],
            [2.0, 0.5, 0.5],
        ])
        modes = np.zeros((1, 4, 3))
        modes[0, 3] = [0.0, 0.0, 0.01]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        result = ana.project_on_internals(0, [[0, 1, 2, 3]])
        assert len(result) == 1
        label, comp, pct = result[0]
        assert label == "D(C0,C1,C2,C3)"
        assert abs(pct - 100.0) < 0.01

    def test_multiple_internals(self):
        atoms = ("C", "C", "C")
        geom = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ])
        modes = np.zeros((1, 3, 3))
        modes[0, 1] = [0.01, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        internals = [[0, 1], [1, 2]]
        result = ana.project_on_internals(0, internals)
        assert len(result) == 2
        label1, comp1, pct1 = result[0]
        label2, comp2, pct2 = result[1]
        assert label1 == "B(C0,C1)"
        assert label2 == "B(C1,C2)"
        assert abs(pct1 + pct2 - 100.0) < 0.01

    def test_zero_displacement_mode(self):
        """Mode with no atomic displacement should give zero projection."""
        atoms = ("C", "O")
        geom = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
        modes = np.zeros((1, 2, 3))
        modes[0, 0] = [1e-6, 0.0, 0.0]
        modes[0, 1] = [1e-6, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        result = ana.project_on_internals(0, [[0, 1]])
        label, comp, pct = result[0]
        assert abs(comp) < 1e-6

    def test_fragment_to_projection_integration(self):
        atoms = ("C", "C", "C")
        geom = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ])
        modes = np.zeros((1, 3, 3))
        modes[0, 0] = [-0.3, 0.0, 0.0]
        modes[0, 1] = [0.3, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        fragments = [[0, 1]]
        internals = NormalModeAnalyzer.derive_internals_from_fragments(fragments)
        result = ana.project_on_internals(0, internals)
        assert len(result) == 1
        label, comp, pct = result[0]
        assert label == "B(C0,C1)"
        assert comp != 0.0
        assert abs(pct - 100.0) < 0.01


class TestDeriveInternalsFromConnectivity:
    def test_diatomic(self):
        symbols = ("C", "O")
        pos = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
        result = NormalModeAnalyzer.derive_internals_from_connectivity(symbols, pos)
        assert result == [[0, 1]]

    def test_triatomic_water(self):
        symbols = ("O", "H", "H")
        pos = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ])
        result = NormalModeAnalyzer.derive_internals_from_connectivity(symbols, pos)
        assert [0, 1] in result
        assert [0, 2] in result
        assert [1, 0, 2] in result
        assert len(result) == 3

    def test_four_atom_chain(self):
        symbols = ("C", "C", "C", "C")
        pos = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.5, 0.0, 0.0],
        ])
        result = NormalModeAnalyzer.derive_internals_from_connectivity(symbols, pos)
        assert [0, 1] in result
        assert [1, 2] in result
        assert [2, 3] in result
        assert [0, 1, 2] in result
        assert [1, 2, 3] in result
        assert [0, 1, 2, 3] in result
        assert len(result) == 6

    def test_nonbonded_atoms(self):
        symbols = ("C", "O")
        pos = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        result = NormalModeAnalyzer.derive_internals_from_connectivity(symbols, pos)
        assert result == []

    def test_single_atom(self):
        result = NormalModeAnalyzer.derive_internals_from_connectivity(("H",), np.zeros((1, 3)))
        assert result == []


class TestProjectOnInternalsWithNormalization:
    def test_percentage_changes_with_full_normalization(self):
        """Percentages should differ when normalising over all internals."""
        atoms = ("O", "C", "O")
        geom = np.array([
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.5, 0.866, 0.0],
        ])
        modes = np.zeros((1, 3, 3))
        modes[0, 2] = [0.0, 0.01, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        # Only the angle
        internals = [[0, 1, 2]]
        pct_only = ana.project_on_internals(0, internals)[0][2]
        assert abs(pct_only - 100.0) < 0.01

        # With full-molecule normalisation (bonds + angle)
        all_internals = NormalModeAnalyzer.derive_internals_from_connectivity(atoms, geom)
        pct_full = ana.project_on_internals(
            0, internals, normalization_internals=all_internals,
        )[0][2]
        assert pct_full < 99.0  # should be less because bonds also contribute

    def test_requested_subset_with_full_normalization(self):
        """Request a subset but normalise over all internals."""
        atoms = ("O", "C", "O")
        geom = np.array([
            [-1.2, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
        ])
        modes = np.zeros((1, 3, 3))
        modes[0, 0] = [-0.3, 0.0, 0.0]
        modes[0, 1] = [0.6, 0.0, 0.0]
        modes[0, 2] = [-0.3, 0.0, 0.0]

        ana = NormalModeAnalyzer(modes, geom, atoms)
        internals = [[0, 1]]  # only one bond
        all_internals = NormalModeAnalyzer.derive_internals_from_connectivity(atoms, geom)

        result = ana.project_on_internals(
            0, internals, normalization_internals=all_internals,
        )
        assert len(result) == 1
        label, comp, pct = result[0]
        assert label == "B(O0,C1)"
        assert pct > 0.0
        assert pct < 100.0  # other internals also contribute


class TestFormatInternalLabel:
    def test_bond_label(self):
        ana = NormalModeAnalyzer(
            np.zeros((1, 2, 3)),
            np.zeros((2, 3)),
            ("C", "O"),
        )
        label = ana._format_internal_label([0, 1])
        assert label == "B(C0,O1)"

    def test_angle_label(self):
        ana = NormalModeAnalyzer(
            np.zeros((1, 3, 3)),
            np.zeros((3, 3)),
            ("O", "C", "O"),
        )
        label = ana._format_internal_label([0, 1, 2])
        assert label == "A(O0,C1,O2)"

    def test_dihedral_label(self):
        ana = NormalModeAnalyzer(
            np.zeros((1, 4, 3)),
            np.zeros((4, 3)),
            ("C", "C", "C", "C"),
        )
        label = ana._format_internal_label([0, 1, 2, 3])
        assert label == "D(C0,C1,C2,C3)"
