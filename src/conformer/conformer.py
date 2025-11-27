
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union
from ase.atoms import Atoms


from src.IOsystem import mkdir
from src.constants import A
from src._calculators.base import BaseCalc

from .energy_data import EnergyRecord, EnergyStore
from .spectral_data import SpectralRecord, SpectralStore

import numpy as np
import random


@dataclass
class Conformer: 
    """
    Dataclass storing all conformer-related data used across the protocol.
    """
    number              : int
    geom                : np.ndarray
    atoms               : tuple
    raw                 : bool            = False

    last_geometry       : np.ndarray      = field(init=False)
    _initial_geometry   : np.ndarray      = field(init=False)
    energies            : EnergyStore     = field(default_factory = EnergyStore)
    active              : bool            = True
    color               : str             = field(default_factory=lambda: "#%06x"%random.randint(0,0xFFFFFF))
    cluster             : Optional[int]   = None
    folder              : str             = field(init=False)
    graphs_data         : SpectralStore   = field(default_factory = SpectralStore)


    def __post_init__(self): 
        self._initial_geometry = self.geom.copy()
        self.last_geometry = self.geom.copy()
        self.folder = f'conf_{self.number}'

        if not self.raw: 
            mkdir(self.folder)

    # ===
    # ASE
    # ===

    def get_ase_atoms(self, calc: BaseCalc) -> Atoms: 
        return Atoms(symbols="".join(tuple(self.atoms)), positions=self.last_geometry, calculator=calc)
    
    # ===
    # Energy helper
    # ===

    def get_energy(self, protocol_number: Union[str, int]):
        energies = self.energies.__getitem__(protocol_number=protocol_number)
        if ~np.isnan(energies.G) or (energies.G is not None):
            return energies.G
        return energies.E

    # ===
    # Properties
    # ===

    @property
    def weight_mass(self):
        return np.sum(
            Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            ).get_masses()
        )

    @property
    def rotatory(self):
        return self.energy.last().B

    @property
    def moment(self):
        return self.energy.last().m

    @property
    def _last_energy(self):
        return self.energies.get_energy()
    
    # ===
    # Geometry helpers
    # ===
    def distance_matrix(self, include_H: bool, geom=None):
        geo = geom if geom is not None else self.last_geometry

        if include_H:
            geo = np.array(geo)
        else:
            mask = self.atoms != "H"
            geo = np.array(geo)[mask]

        return np.linalg.norm(geo[:, None, :] - geo[None, :, :], axis=-1)

    # ===
    # Deserialization
    # === 

    @staticmethod
    def load_raw(data):
        c = Conformer(
            number=data["number"],
            geom=data["last_geometry"],
            atoms=data["atoms"],
            raw=True,
        )
        c.energies = data["energies"]
        c.active = data["active"]
        return c
    
    # === 
    # Sorting support
    # ===
    def __lt__(self, other):
        if not self.active:
            return 0 < other._last_energy
        return self._last_energy < other._last_energy

    def __gt__(self, other):
        if not self.active:
            return 0 > other._last_energy
        return self._last_energy > other._last_energy

    def __eq__(self, other):
        if not self.active:
            return 0 == other._last_energy
        return self._last_energy == other._last_energy


if __name__=='__main__':
    Conformer() 