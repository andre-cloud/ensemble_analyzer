
from typing import List
from pathlib import Path
import json
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer.io_utils import SerialiseEncoder



class ProtocolManager:
    """
    Manages protocol persistence and loading.
    
    Responsibilities:
    - Protocol serialization
    - Protocol validation
    - Last completed protocol tracking
    """
    
    def __init__(self, protocol_file: str = "protocol_dump.json", last_protocol_file: str = "last_protocol"):
        
        self.protocol_file = Path(protocol_file)
        self.last_protocol_file = Path(last_protocol_file)
    
    def save(self, protocols: List[Protocol]) -> None:

        data = {p.number: p.__dict__ for p in protocols}
        with open(self.protocol_file, 'w') as f:
            json.dump(data, f, indent=4, cls=SerialiseEncoder)
    
    def load(self) -> List[Protocol]:

        with open(self.protocol_file) as f:
            data = json.load(f)
        return [Protocol(**data[key]) for key in data]
    
    def save_last_completed(self, protocol_number: int) -> None:

        with open(self.last_protocol_file, 'w') as f:
            f.write(str(protocol_number))
    
    def load_last_completed(self) -> int:

        with open(self.last_protocol_file) as f:
            return int(f.read().strip())