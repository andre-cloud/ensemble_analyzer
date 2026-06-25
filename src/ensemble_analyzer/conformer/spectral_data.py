from dataclasses import dataclass, field
from typing import Optional, Dict, Union, Literal

import numpy as np

from collections import defaultdict


@dataclass
class SpectralRecord:
    """
    Data container for spectral transitions (impulses).
    """

    X : np.ndarray # energy impulses
    Y : np.ndarray # impulse intensity

    def __post_init__(self) -> None:
        """Validate and convert X/Y to 1D numpy arrays of matching shape."""
        if not isinstance(self.X, np.ndarray):
            self.X = np.array(self.X)
        if not isinstance(self.Y, np.ndarray):
            self.Y = np.array(self.Y)

        if self.X.shape != self.Y.shape:
            raise ValueError(
                f"X and Y must have same shape. Got X: {self.X.shape}, Y: {self.Y.shape}"
            )
        
        if self.X.ndim != 1:
            raise ValueError(f"X and Y must be 1D arrays. Got {self.X.ndim}D")

    def as_dict(self) -> dict:
        """Convert to serializable dictionary."""

        return {
            "X": self.X.tolist(),
            "Y": self.Y.tolist(),
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'SpectralRecord':
        """Reconstruct from dictionary."""

        return cls(
            X=np.array(data["X"]),
            Y=np.array(data["Y"])
        )
    
    def __len__(self) -> int:
        """Number of spectral transitions."""
        return len(self.X)
    
    @property
    def is_empty(self) -> bool:
        """Check whether the record contains no transitions."""
        return len(self.X) == 0

    

@dataclass
class SpectralStore:
    """
    Hierarchical storage for spectral data: Protocol -> GraphType -> SpectralRecord.
    Example: store[1]['IR'] -> SpectralRecord(...)
    """

    data: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(SpectralRecord)))


    def add(self, protocol_number: int, graph_type: Literal['IR', 'VCD', 'UV', 'ECD'], record: SpectralRecord) -> None:
        """Store a spectral record for the given protocol and graph type."""
        self.data[int(protocol_number)][str(graph_type)] = record

    def __getitem__(self, key) -> SpectralRecord:
        """Retrieve a spectral record by (protocol_number, graph_type) tuple."""
        protocol_number, graph_type = key
        return self.data[int(protocol_number)][str(graph_type)]

    def __contains__(self, protocol_number: int) -> bool:
        """Check if data exists for the given protocol number."""
        return int(protocol_number) in self.data
    
    def has_graph_type(self, protocol_number: int, graph_type: Literal['IR', 'VCD', 'UV', 'ECD']) -> bool:
        """Check if a specific graph type exists for the given protocol."""
        return graph_type in self.data[int(protocol_number)]

    def as_dict(self) -> dict:
        """Used for checkpoint serialization."""
        return {k: {k1: v1.as_dict() for k1, v1 in v.items()} for k, v in self.data.items()}
    
    def load(self, input_dict: dict) -> None:
        """Restore the store from a serialized dictionary."""
        self.data = defaultdict(lambda: defaultdict(SpectralRecord))
        for proto_str, graphs in input_dict.items():
            proto = int(proto_str)
            for graph_type, record_dict in graphs.items():
                self.data[proto][graph_type] = SpectralRecord.from_dict(record_dict)