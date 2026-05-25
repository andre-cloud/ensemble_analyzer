
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union
from ase.atoms import Atoms


from ensemble_analyzer.io_utils import mkdir
from ensemble_analyzer.calculators.base import BaseCalc

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


    def __post_init__(self) -> None:
        """Initialize derived fields from the input geometry."""
        self._initial_geometry = self.geom.copy()
        self.last_geometry = self.geom.copy()
        self.folder = f'conf_{self.number}'

        if not self.raw: 
            mkdir(self.folder)

    # ===
    # ASE
    # ===

    def get_ase_atoms(self, calc: BaseCalc) -> Atoms:
        """Build an ASE Atoms object with the given calculator."""
        return Atoms(symbols="".join(tuple(self.atoms)), positions=self.last_geometry, calculator=calc)
    
    # ===
    # Energy helper
    # ===

    def get_energy(self, protocol_number: int) -> float:
        """Return Gibbs free energy if available, else electronic energy."""
        energies = self.energies[protocol_number]
        if not np.isnan(energies.G):
            return energies.G
        return energies.E
    
    def create_log(self, protocol_number: int, monitor_internals: List[List[int]]) -> tuple:
        """Build a log tuple with conformer data and optional internal coordinates."""
        e, g_e, g, b, erel, pop, time = self.energies.log_info(protocol_number=protocol_number)

        monitor : List[float] = []
        if len(monitor_internals) > 0:
            atoms = Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            )
            for internal in monitor_internals:
                if len(internal) == 2:
                    monitor.append(float(atoms.get_distance(*internal)))
                if len(internal) == 3:
                    monitor.append(float(atoms.get_angle(*internal)))
                if len(internal) == 4:
                    monitor.append(float(atoms.get_dihedral(*internal)))
                    
        if len(monitor)==0:
            return self.number, e, g_e, g, b, erel, pop, time, self.cluster
        
        return  self.number, e, g_e, g, b, erel, pop, time, self.cluster, *monitor


    
    def write_xyz(self) -> str:
        """Generate XYZ-formatted string for file output.

        Returns:
            str: XYZ-formatted string, or empty string if conformer is inactive.
        """
        
        if not self.active:
            return ""

        # Header
        header = f'{len(self.atoms)}\nCONFORMER {self.number} {self._last_energy:10f}'

        # Atoms and positions
        atom_lines = [
            f"{a}  {x:14.6f}  {y:14.6f}  {z:14.6f}"
            for a, (x, y, z) in zip(self.atoms, self.last_geometry)
        ]

        txt = "\n".join([header] + atom_lines)

        return txt

    # ===
    # Properties
    # ===

    @property
    def weight_mass(self) -> float:
        """Total mass of the conformer in atomic mass units."""
        return np.sum(
            Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            ).get_masses()
        )

    @property
    def rotatory(self) -> float:
        """Rotational constant norm from the most recent protocol step."""
        return self.energies.last().B

    @property
    def moment(self) -> float:
        """Dipole moment norm from the most recent protocol step."""
        return self.energies.last().m

    @property
    def _last_energy(self) -> float:
        """Gibbs or electronic energy from the most recent protocol step."""
        return self.energies.get_energy()
    
    # ===
    # Geometry helpers
    # ===
    def distance_matrix(self, include_H: bool, geom: Optional[np.ndarray] = None) -> np.ndarray:
        """Compute the pairwise distance matrix for the conformer.

        Args:
            include_H: Whether to include hydrogen atoms.
            geom: Optional geometry array. Uses last_geometry if None.

        Returns:
            np.ndarray: Pairwise distance matrix.
        """
        geo = geom if geom is not None else self.last_geometry

        if include_H:
            geo = np.array(geo)
        else:
            mask = np.array(self.atoms) != "H"
            geo = np.array(geo)[mask]

        return np.linalg.norm(geo[:, None, :] - geo[None, :, :], axis=-1)

    # ===
    # Deserialization
    # === 

    @staticmethod
    def load_raw(data: dict) -> 'Conformer':
        """Deserialize a Conformer from a dictionary.

        Args:
            data: Dictionary with conformer data (number, last_geometry, atoms,
                  energies, graphs_data, active).

        Returns:
            Conformer: Restored conformer instance.
        """
        c = Conformer(
            number=data["number"],
            geom=data["last_geometry"],
            atoms=data["atoms"],
            raw=True,
        )
        
        c.energies.load(data["energies"])
        c.graphs_data.load(data["graphs_data"])
        
        c.active = data["active"]
        return c
    
    # === 
    # Sorting support
    # ===
    def __lt__(self, other: 'Conformer') -> bool:
        """Compare conformers by energy for sorting (lowest first)."""
        if not self.active:
            return 0 < other._last_energy
        return self._last_energy < other._last_energy

    def __gt__(self, other: 'Conformer') -> bool:
        """Compare conformers by energy for sorting (highest first)."""
        if not self.active:
            return 0 > other._last_energy
        return self._last_energy > other._last_energy

    def __eq__(self, other: 'Conformer') -> bool:
        """Check if two conformers have the same energy."""
        if not self.active:
            return 0 == other._last_energy
        return self._last_energy == other._last_energy