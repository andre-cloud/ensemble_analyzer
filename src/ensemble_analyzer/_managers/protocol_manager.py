
from typing import List
from pathlib import Path
import json
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer.io_utils import SerialiseEncoder



class ProtocolManager:
    """
    Manages protocol persistence and loading.
    
    Responsibilities:
    - Protocol serialization to JSON.
    - Protocol validation/deserialization.
    - Tracking the last completed protocol step for restart capability.
    """
    
    def __init__(self, protocol_file: str = "protocol_dump.json", last_protocol_file: str = "last_protocol"):
        """
        Initialize the manager.

        Args:
            protocol_file (str): Filename for the dumped protocol JSON.
            last_protocol_file (str): Filename tracking the last finished step.
        """

        self.protocol_file = Path(protocol_file)
        self.last_protocol_file = Path(last_protocol_file)
    
    def save(self, protocols: List[Protocol]) -> None:
        """
        Serialize and save the list of protocols to disk.

        Used to persist the exact protocol configuration used for a run.

        Args:
            protocols (List[Protocol]): List of Protocol objects to save.

        Returns:
            None
        """

        data = {p.number: p.__dict__ for p in protocols}
        with open(self.protocol_file, 'w') as f:
            json.dump(data, f, indent=4, cls=SerialiseEncoder)
    
    def load(self) -> List[Protocol]:
        """
        Load protocols from the dump file.

        Returns:
            List[Protocol]: Reconstructed list of Protocol objects.
        """

        with open(self.protocol_file) as f:
            data = json.load(f)
        return [Protocol(**data[key]) for key in data]
    
    def save_last_completed(self, protocol_number: int) -> None:
        """
        Update the marker for the last successfully completed protocol.

        Args:
            protocol_number (int): The ID of the completed protocol.

        Returns:
            None
        """

        with open(self.last_protocol_file, 'w') as f:
            f.write(str(protocol_number))
    
    def load_last_completed(self) -> int:
        """
        Retrieve the ID of the last completed protocol.

        Returns:
            int: Protocol ID (returns 0 if file doesn't exist).
        """

        with open(self.last_protocol_file) as f:
            return int(f.read().strip())