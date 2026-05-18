"""
Tests for Checkpoint Manager.
Verifies saving and loading of calculation state.
"""

import pytest
import json
from unittest.mock import MagicMock, patch, mock_open
from ensemble_analyzer._managers.checkpoint_manager import CheckpointManager

class TestCheckpointManager:
    
    @pytest.fixture
    def manager(self):
        # FIX: Argument name is checkpoint_file, not filename
        return CheckpointManager(checkpoint_file="checkpoint.json")

    def test_save_checkpoint(self, manager, mock_logger, mock_conformer):
        """Test saving conformers to JSON."""
        mock_conformer.number = 1
        mock_conformer.active = True
        
        # We need to ensure the parent directory exists logic doesn't crash the mock
        with patch("ensemble_analyzer._managers.checkpoint_manager.tempfile.NamedTemporaryFile") as mock_temp:
            # Setup temp file mock context manager
            mock_temp_file = MagicMock()
            mock_temp.return_value.__enter__.return_value = mock_temp_file
            mock_temp_file.name = "temp_checkpoint.tmp"
            
            with patch("ensemble_analyzer._managers.checkpoint_manager.shutil.move") as mock_move:
                with patch("json.dump") as mock_json_dump:
                    
                    manager.save([mock_conformer], mock_logger, log=True)
                    
                    # Verify json dump called with correct data structure
                    mock_json_dump.assert_called()
                    args, _ = mock_json_dump.call_args
                    data_arg = args[0]
                    assert 1 in data_arg # Key should be conformer number
                    
                    # Verify move
                    mock_move.assert_called_with("temp_checkpoint.tmp", "checkpoint.json")
                    
                    # Verify logging
                    mock_logger.checkpoint_saved.assert_called()

    def test_load_checkpoint_success(self, manager):
        """Test loading a valid checkpoint."""
        fake_data = {
            "1": {
                "number": 1,
                "atoms": ["C"],
                "last_geometry": [[0,0,0]],
                "active": True,
                "energies": {},
                "graphs_data": {}
            }
        }
        
        # CheckpointManager.load uses Path.exists(), then open().
        with patch("pathlib.Path.exists", return_value=True):
            with patch("builtins.open", mock_open(read_data=json.dumps(fake_data))):
                with patch("ensemble_analyzer._managers.checkpoint_manager.Conformer") as MockConfClass:
                    MockConfClass.load_raw.return_value = MagicMock(number=1)
                    
                    confs = manager.load()
                    
                    assert len(confs) == 1
                    MockConfClass.load_raw.assert_called()

    def test_load_checkpoint_not_found(self, manager):
        """Test loading when file doesn't exist."""
        with patch("pathlib.Path.exists", return_value=False):
            with pytest.raises(FileNotFoundError):
                manager.load()