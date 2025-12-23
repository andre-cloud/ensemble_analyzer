"""
Tests for the Ensemble Input/Output module.
Covers parsing of XYZ strings, reading ensemble files, and saving snapshots.
"""

import pytest
import numpy as np
from unittest.mock import mock_open, patch
from ensemble_analyzer.ensemble_io import read_ensemble, save_snapshot, _parse_xyz_str
from unittest.mock import MagicMock

class TestEnsembleIO:
    
    def test_parse_xyz_str_standard(self):
        """Test parsing of a standard XYZ block without raw energy."""
        lines = ["2","","C 0.0 0.0 0.0", "O 0.0 0.0 1.2"]
        atoms, geom, energy = _parse_xyz_str(lines, raw=False)
        assert np.array_equal(atoms, np.array(["C", "O"]))
        assert geom.shape == (2, 3)
        assert energy is None

    def test_parse_xyz_str_raw_energy(self):
        """Test parsing of raw energy from the comment line."""
        lines = ["3", "Comment line with energy -150.12345", "C 0.0 0.0 0.0"]
        atoms, geom, energy = _parse_xyz_str(lines, raw=True)
        assert energy == -150.12345

    def test_read_ensemble_valid(self, mock_logger):
        """Test reading a valid XYZ file containing multiple conformers."""
        content = "3\nTest\nC 0 0 0\nO 0 0 1\nH 0 1 0\n3\nTest2\nC 1 0 0\nO 1 0 1\nH 1 1 0"
        with patch("builtins.open", mock_open(read_data=content)):
            with patch("ensemble_analyzer.ensemble_io.Conformer") as MockConf:
                res = read_ensemble("file.xyz", mock_logger)
                assert len(res) == 2
                assert MockConf.call_count == 2

    def test_read_ensemble_invalid_extension(self, mock_logger):
        """Test that reading a non-XYZ file raises a ValueError."""
        with pytest.raises(ValueError) as excinfo:
            read_ensemble("test.txt", mock_logger)
        assert "must be an XYZ" in str(excinfo.value)

    def test_save_snapshot_filtering(self, mock_logger, mock_conformer):
        """Test that save_snapshot only writes valid/active conformers."""
        bad_conf = MagicMock()
        bad_conf.write_xyz.return_value = None
        
        confs = [mock_conformer, bad_conf]
        
        with patch("builtins.open", mock_open()) as mocked_file:
            save_snapshot("out.xyz", confs, mock_logger)
            handle = mocked_file()
            # Ensure only the valid conformer's data is written
            handle.write.assert_called_once_with(mock_conformer.write_xyz.return_value)