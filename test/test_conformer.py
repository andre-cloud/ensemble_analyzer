"""
Tests for the Conformer dataclass.
Verifies initialization, XYZ generation, sorting, and distance matrix calculation.
"""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from ensemble_analyzer._conformer.conformer import Conformer

class TestConformer:
    
    @pytest.fixture
    def basic_conformer(self):
        """Fixture for a basic Conformer, mocking filesystem operations."""
        with patch("ensemble_analyzer._conformer.conformer.mkdir"):
            return Conformer(
                number=1,
                geom=np.array([[0.,0.,0.], [0.,0.,1.]]),
                atoms=("C", "O")
            )

    def test_initialization(self, basic_conformer):
        """Test correct attribute initialization."""
        assert basic_conformer.folder == "conf_1"
        assert basic_conformer.active is True
        assert basic_conformer.atoms == ("C", "O")

    def test_write_xyz(self, basic_conformer):
        """Test generation of XYZ string format."""
        basic_conformer.energies.get_energy = MagicMock(return_value=-150.0)
        
        xyz = basic_conformer.write_xyz()
        assert "2" in xyz.splitlines()[0] # Number of atoms
        assert "CONFORMER 1" in xyz
        assert "C" in xyz and "O" in xyz

    def test_inactive_xyz_empty(self, basic_conformer):
        """Test that inactive conformers return an empty XYZ string."""
        basic_conformer.active = False
        assert basic_conformer.write_xyz() == ""

    def test_sorting(self, basic_conformer):
        """Test conformer comparison based on energy."""
        c1 = basic_conformer
        c1.energies.get_energy = MagicMock(return_value=-100.0)
        
        c2 = Conformer(2, np.array([]), ("H",), raw=True)
        c2.energies.get_energy = MagicMock(return_value=-200.0)
        
        # Lower energy means "less than" in sorting
        assert c2 < c1
        assert c1 > c2

    def test_distance_matrix(self, basic_conformer):
        """Test calculation of the interatomic distance matrix."""
        dm = basic_conformer.distance_matrix(include_H=True)
        assert dm[0,1] == 1.0 # Distance between (0,0,0) and (0,0,1)
        assert dm[1,0] == 1.0
        assert dm[0,0] == 0.0