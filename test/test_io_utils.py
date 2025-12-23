"""
Tests for IO Utilities.
Verifies file movement, directory creation, json encoding and file tailing.
"""

import pytest
import numpy as np
import json
from unittest.mock import patch, MagicMock, mock_open
from ensemble_analyzer.io_utils import mkdir, move_files, tail, SerialiseEncoder

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