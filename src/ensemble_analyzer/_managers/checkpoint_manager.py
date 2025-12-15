
from typing import List
from pathlib import Path
import json
import tempfile
import shutil

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._logger.logger import Logger
from ensemble_analyzer.io_utils import SerialiseEncoder

class CheckpointManager:
    """
    Manages atomic checkpoint saves and loads.
    
    Ensures data integrity by using atomic file operations (write to temp -> move)
    to prevent corruption if the program is interrupted during a save.
    """
    
    def __init__(self, checkpoint_file: str = "checkpoint.json"):
        """
        Initialize the manager.

        Args:
            checkpoint_file (str): Path to the checkpoint JSON file.
        """

        self.checkpoint_file = Path(checkpoint_file)
    
    def save(self, ensemble: List[Conformer], logger: Logger, log: bool = False) -> None:
        """
        Save the ensemble state atomically.

        Args:
            ensemble (List[Conformer]): The list of conformers to persist.
            logger (Logger): Logger instance for tracking.
            log (bool): Whether to write a log message confirming the save.

        Returns:
            None
        """
        
        data = {conf.number: conf.__dict__ for conf in ensemble}
        
        # Write to temporary file first
        with tempfile.NamedTemporaryFile(
            'w',
            delete=False,
            suffix='.tmp',
            dir=self.checkpoint_file.parent
        ) as tmp:
            json.dump(data, tmp, indent=4, cls=SerialiseEncoder)
            tmp_path = Path(tmp.name)
        
        # Atomic move
        shutil.move(str(tmp_path), str(self.checkpoint_file))
        
        if log:
            logger.checkpoint_saved(conformer_count=len(ensemble))
    
    def load(self) -> List[Conformer]:
        """
        Load the ensemble from the checkpoint file.

        Returns:
            List[Conformer]: The restored ensemble.
            
        Raises:
            FileNotFoundError: If the checkpoint file does not exist.
        """
        
        if not self.checkpoint_file.exists():
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_file}")
        
        with open(self.checkpoint_file) as f:
            data = json.load(f)
        
        return [Conformer.load_raw(data[key]) for key in data]