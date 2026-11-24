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

    def test_save_checkpoint(self, manager, mock_logger, mock_conformer, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        """Test saving conformers to JSON."""
        mock_conformer.number = 1
        mock_conformer.active = True
        
        with patch("ensemble_analyzer._managers.checkpoint_manager.write_json",
                   side_effect=lambda data, f, indent: f.write("{}")) as mock_write_json:
            manager.save([mock_conformer], mock_logger, log=True)

            # Verify write_json called with correct data structure
            args, _ = mock_write_json.call_args
            assert '1' in args[0]  # Key should be conformer number

            # Atomic replace: the checkpoint exists, the temporary file does not
            assert manager.checkpoint_file.exists()
            assert not manager.checkpoint_file.with_suffix(".tmp").exists()

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