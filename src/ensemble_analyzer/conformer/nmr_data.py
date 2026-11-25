"""Per-conformer NMR results: isotropic shieldings and spin-spin couplings, per protocol step.

Atom indices are 0-based; couplings are stored once per pair (i < j), in Hz, for the isotopes
requested in the step.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple


@dataclass
class NMRRecord:
    shielding: Dict[int, float] = field(default_factory=dict)   # sigma_iso, ppm
    aniso: Dict[int, float] = field(default_factory=dict)       # anisotropy, ppm
    coupling: Dict[Tuple[int, int], float] = field(default_factory=dict)  # J, Hz

    def __post_init__(self):
        self.shielding = {int(k): float(v) for k, v in self.shielding.items()}
        self.aniso = {int(k): float(v) for k, v in self.aniso.items()}
        self.coupling = {(min(int(a), int(b)), max(int(a), int(b))): float(v)
                         for (a, b), v in self.coupling.items()}

    def j(self, a: int, b: int) -> Optional[float]:
        return self.coupling.get((min(a, b), max(a, b)))

    def has_data(self) -> bool:
        return bool(self.shielding or self.coupling)

    def as_dict(self) -> dict:
        return {
            "shielding": {str(k): v for k, v in self.shielding.items()},
            "aniso": {str(k): v for k, v in self.aniso.items()},
            "coupling": [[a, b, v] for (a, b), v in sorted(self.coupling.items())],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "NMRRecord":
        return cls(
            shielding=data.get("shielding", {}),
            aniso=data.get("aniso", {}),
            coupling={(int(a), int(b)): float(v) for a, b, v in data.get("coupling", [])},
        )


@dataclass
class NMRStore:
    data: Dict[int, NMRRecord] = field(default_factory=dict)

    def add(self, protocol_number: int, record: NMRRecord) -> None:
        self.data[int(protocol_number)] = record

    def get(self, protocol_number) -> Optional[NMRRecord]:
        return self.data.get(int(protocol_number))

    def __contains__(self, protocol_number) -> bool:
        return int(protocol_number) in self.data

    def as_dict(self) -> dict:
        return {str(k): v.as_dict() for k, v in self.data.items()}

    def load(self, input_dict: dict) -> None:
        self.data = {int(k): NMRRecord.from_dict(v) for k, v in (input_dict or {}).items()}
