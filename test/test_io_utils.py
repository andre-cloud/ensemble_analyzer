"""
Tests for IO Utilities.
Verifies file movement, directory creation, json encoding and file tailing.
"""

import io
import json

import numpy as np
import pytest
from unittest.mock import patch, MagicMock, mock_open

from ensemble_analyzer.io_utils import mkdir, move_files, tail, write_json, SerialiseEncoder, _serialise
from ensemble_analyzer.protocol.solvent import Solvent

class TestIOUtils:

    @patch("ensemble_analyzer.io_utils.Path")
    def test_mkdir(self, mock_path):
        """Test directory creation logic."""
        # Mock Path object return
        mock_dir = MagicMock()
        mock_path.return_value = mock_dir
        
        # Execute
        res = mkdir("test_dir")
        
        # Assert
        assert res is True
        mock_path.assert_called_with("test_dir")
        mock_dir.mkdir.assert_called_with(parents=True, exist_ok=True)

    @patch("ensemble_analyzer.io_utils.shutil")
    @patch("ensemble_analyzer.io_utils.Path")
    def test_move_files(self, mock_path, mock_shutil, mock_conformer, mock_protocol):
        """Test moving files based on label."""
        mock_conformer.folder = "conf_1"
        mock_conformer.number = 1
        mock_protocol.number = 1
        label = "opt"
        
        # Mock current working directory and file listing
        mock_cwd = MagicMock()
        mock_path.cwd.return_value = mock_cwd
        
        # Create a mock file that matches the label
        mock_file = MagicMock()
        mock_file.name = "opt_output.log"
        # Mock iterator to return our file
        mock_cwd.iterdir.return_value = [mock_file]
        
        # Execute
        move_files(mock_conformer, mock_protocol, label)
        
        # Assert move was called
        # dest path construction involves joins, checking exact string is complex with mocks
        # checking called is sufficient for logic flow
        assert mock_shutil.move.called

    def test_tail(self):
        """Test reading the last N lines of a file."""
        content = "Line 1\nLine 2\nLine 3\nLine 4"
        
        with patch("pathlib.Path.open", mock_open(read_data=content)):
            # Read last 2 lines
            result = tail("dummy.log", num_lines=2)
            assert result == "Line 3\nLine 4"

    def test_serialise_encoder(self):
        """Test custom JSON encoder for NumPy arrays and Objects."""
        
        # Test NumPy array serialization
        data = {"array": np.array([1, 2, 3])}
        json_str = json.dumps(data, cls=SerialiseEncoder)
        assert "[1, 2, 3]" in json_str
        
        # Test Object serialization (via __dict__)
        class DummyObj:
            def __init__(self):
                self.val = 42
                
        data_obj = {"obj": DummyObj()}
        json_str_obj = json.dumps(data_obj, cls=SerialiseEncoder)
        assert '"val": 42' in json_str_obj

    def test_serialise_int_keys(self):
        """_serialise converts integer dict keys to strings."""
        data = {0: {"E": -1961.03, "G": -1960.48}, 1: {"E": -1960.12}}
        result = _serialise(data)
        assert "0" in result and "1" in result
        assert all(isinstance(k, str) for k in result)

    def test_serialise_np_int_keys(self):
        """_serialise converts numpy integer dict keys to strings."""
        data = {np.int64(0): {"E": -1961.03}}
        result = _serialise(data)
        assert "0" in result
        assert all(isinstance(k, str) for k in result)

    def test_write_json_int_keys(self):
        """write_json produces valid JSON with integer keys."""
        data = {0: {"E": -1961.03}, 2: {"E": -1960.48}}
        buf = io.StringIO()
        write_json(data, buf)
        result = json.loads(buf.getvalue())
        assert result["0"]["E"] == -1961.03
        assert result["2"]["E"] == -1960.48

    def test_write_json_dataclass(self):
        """write_json serialises dataclass objects like Solvent."""
        data = {
            "0": {
                "functional": "r2scan-3c",
                "solvent": Solvent("CHCl3"),
            }
        }
        buf = io.StringIO()
        write_json(data, buf)
        result = json.loads(buf.getvalue())
        assert result["0"]["solvent"] == {"solvent": "CHCl3", "smd": False}
        assert result["0"]["functional"] == "r2scan-3c"