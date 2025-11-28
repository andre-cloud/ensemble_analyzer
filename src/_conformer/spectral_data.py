from dataclasses import dataclass, field
from typing import Optional, Dict, Union, Literal

import numpy as np

from collections import defaultdict


@dataclass
class SpectralRecord:
    """
    Storing the impulses
    """
    X : np.ndarray # energy impulses
    Y : np.ndarray # impulse intensity

    def as_dict(self):
        return {"x": self.X, "y": self.Y}
    

@dataclass
class SpectralStore:
    data: Dict = field(default_factory=lambda: defaultdict(lambda: defaultdict(SpectralRecord)))


    def add(self, protocol_number:int, graph_type: Literal['IR', 'VCD', 'UV', 'ECD'], record: SpectralRecord):
        self.data[int(protocol_number)][str(graph_type)] = record

    def __getitem__(self, protocol_number:int, graph_type: Union[int, str]) -> SpectralRecord:
        return self.data[int(protocol_number)][str(graph_type)]

    def __contains__(self, protocol_number:int) -> bool:
        return int(protocol_number) in self.data
    
    def __has_graph_type__(self, protocol_number:int, graph_type: Literal['IR', 'VCD', 'UV', 'ECD']):
        return graph_type in self.data[int(protocol_number)]

    def as_dict(self):
        """Used for checkpoint serialization"""
        return {k: {k1: v1} for k, v in self.data.items() for k1, v1 in v.items()}
    
    def load(self, input_dict: Dict[int, Dict[str, SpectralRecord]]):
        self.data = defaultdict(lambda: defaultdict(SpectralRecord))  # reset self.data
        for proto_str, graphs in input_dict.get('data', {}).items():
            proto = int(proto_str)
            for graph_type, record_dict in graphs.items():
                x_array = np.array(record_dict['X'])
                y_array = np.array(record_dict['Y'])
                self.data[proto][graph_type] = SpectralRecord(X=x_array, Y=y_array)