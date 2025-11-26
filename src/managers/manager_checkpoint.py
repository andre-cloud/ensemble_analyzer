
from typing import List
from pathlib import Path
import json

from src.conformer import Conformer
from src.logger.logger import Logger
from src.IOsystem import SerialiseEncoder

class CheckpointManager:
    """
    Manages atomic checkpoint saves and loads.
    
    Responsibilities:
    - Atomic file writes to prevent corruption
    - JSON serialization with proper encoding
    - Conformer state persistence
    """
    
    def __init__(self, checkpoint_file: str = "checkpoint.json"):
        self.checkpoint_file = Path(checkpoint_file)
    
    def save(self, ensemble: List[Conformer], logger: Logger) -> None:
        """
        Save checkpoint atomically.
        
        Args:
            ensemble: List of conformers to save
            logger: Logger for recording operation
        """
        import tempfile
        import shutil
        
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
        
        # Log checkpoint
        file_size_mb = self.checkpoint_file.stat().st_size / (1024 * 1024)
        logger.checkpoint_saved(
            conformer_count=len(ensemble),
            file_size_mb=file_size_mb
        )
    
    def load(self) -> List[Conformer]:
        """
        Load checkpoint from file.
        
        Returns:
            List of conformers
            
        Raises:
            FileNotFoundError: If checkpoint doesn't exist
        """
        if not self.checkpoint_file.exists():
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_file}")
        
        with open(self.checkpoint_file) as f:
            data = json.load(f)
        
        return [Conformer.load_raw(data[key]) for key in data]