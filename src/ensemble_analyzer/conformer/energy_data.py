from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, Tuple, Union, TYPE_CHECKING
import numpy as np
from ensemble_analyzer.constants import ROT_CONST_FACTOR
from ensemble_analyzer.rrho import free_gibbs_energy
from ase import Atoms

if TYPE_CHECKING:
    from ensemble_analyzer.conformer.conformer import Conformer



@dataclass
class EnergyRecord:
    """
    Data container for energetic and thermodynamic properties of a conformer.
    """

    E       : float                     = 0.0           # Electronic Energy [Eh]
    G       : float                     = np.nan        # Gibbs Free Energy [Eh]
    H       : float                     = np.nan        # Enthalpy [Eh]
    S       : float                     = np.nan        # Total Entropy [Eh]
    G_E     : float                     = np.nan        # Thermal Correction to Gibbs (G-E)
    zpve    : float                     = np.nan        # Zero Point Vibrational Energy
    B       : Optional[float]           = None          # Rotational Constant Norm [cm-1]
    B_vec   : Optional[np.ndarray]      = None          # Rotational Constants Vector [cm-1]
    m       : Optional[float]           = None          # Dipole Moment Norm [Debye]
    m_vec   : Optional[np.ndarray]      = None          # Dipole Moment Vector
    Pop     : float                     = np.nan        # Boltzmann Population [%]
    time    : Optional[float]           = None          # Calculation elapsed time [s]
    Erel    : float                     = np.nan        # Relative Energy [kcal/mol]
    Freq    : Optional[np.ndarray]      = None          # Vibrational Frequencies [cm-1]
    NormalModes : Optional[np.ndarray]  = None          # Normal mode vectors (n_modes, n_atoms, 3)
    calculator : str                     = ""            # Calculator name (e.g. "orca", "tblite")

    def as_dict(self) -> dict:
        """Convert record to dictionary, handling numpy arrays."""

        data = asdict(self)

        if data['Freq'] is not None:
            data['Freq'] = np.array(data['Freq'])
            neg_mask = data['Freq'] < 0
            if data['NormalModes'] is not None:
                normal_modes = np.array(data['NormalModes'])
                if len(normal_modes) == len(neg_mask):
                    normal_modes = normal_modes[neg_mask]
                data['NormalModes'] = normal_modes.tolist()
            data['Freq'] = data['Freq'].tolist()

        for key in ['B_vec', 'm_vec']:
            if data[key] is not None:
                data[key] = data[key].tolist()
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> 'EnergyRecord':
        """Create record from dictionary, restoring numpy arrays."""

        for key in ['B_vec', 'm_vec', 'Freq', 'NormalModes']:
            if data.get(key) is not None:
                data[key] = np.array(data[key])
        
        return cls(**data)


@dataclass
class EnergyStore:
    """
    Dictionary-like store for EnergyRecords indexed by protocol number.
    """
    
    data: Dict[str, EnergyRecord] = field(default_factory=dict)

    def add(self, protocol_number: str, record: EnergyRecord) -> None:
        """Add a record for a specific protocol step."""

        self.data[str(protocol_number)] = record

    def last(self) -> EnergyRecord:
        """Retrieve the record from the most recent protocol step."""
        
        if not self.data:
            return EnergyRecord()
        last_key = list(self.data.keys())[-1]
        return self.data[last_key]

    def __getitem__(self, protocol_number: str) -> 'EnergyRecord':
        """Retrieve the record for a given protocol number, or an empty record."""
        if self.__contains__(protocol_number=protocol_number):
            return self.data.get(str(protocol_number))
        
        return EnergyRecord()

    def __contains__(self, protocol_number: str) -> bool:
        """Check if a record exists for the given protocol number."""
        return str(protocol_number) in self.data

    def as_dict(self) -> dict:
        """Serialize to a dictionary for checkpoint storage."""
        return {k: v.as_dict() for k, v in self.data.items()}
    
    def get_energy(self, protocol_number: str = None) -> float:
        """Return Gibbs free energy from a protocol, else electronic energy."""
        data = self[protocol_number] if protocol_number is not None else self.last()
        if not np.isnan(data.G):
            return data.G
        return data.E
    
    def set(self, protocol_number: str, property: str, value: Union[float, np.ndarray]) -> None:
        """Set a specific property on an existing EnergyRecord.

        Args:
            protocol_number: Protocol step number.
            property: Attribute name on EnergyRecord (e.g. 'E', 'G', 'Pop').
            value: Value to set.

        Raises:
            KeyError: If no record exists for the given protocol.
            AttributeError: If the property does not exist on EnergyRecord.
        """
        protocol_number = str(protocol_number)
        if not self.__contains__(protocol_number):
            raise KeyError(f"Protocol {protocol_number} not found in EnergyStore")
        
        if not hasattr(self.data[protocol_number], property):
            raise AttributeError(
                f"EnergyRecord has no attribute '{property}'. "
                f"Valid: E, G, H, S, G_E, zpve, B, B_vec, m, m_vec, Pop, time, Erel, Freq, calculator"
            )
        
        setattr(self.data[protocol_number], property, value)
    
    def log_info(self, protocol_number: str) -> Tuple[float]:
        """Format energy data for log output."""
        data = self.__getitem__(str(protocol_number))
        erel = f'{data.Erel:.2f}' if not np.isnan(data.Erel) else np.nan
        pop = f'{data.Pop:.2f}' if not np.isnan(data.Pop) else np.nan

        b_str = f'{data.B:.5f}' if data.B is not None else 'N/A'
        time_str = f'{data.time:.2f}' if data.time is not None else 'N/A'
        return data.E, data.G_E, data.G, b_str, erel, pop, time_str

    def load(self, input_dict: dict) -> None:
        """Restore the store from a serialized dictionary."""
        self.data = dict()
        for proto_str, vals in input_dict.get('data', {}).items():
            proto = str(proto_str)
                        
            self.data[proto] = EnergyRecord.from_dict(data=vals)

    def get_last_freq(self, protocol_number: str) -> np.ndarray:
        """Retrieve frequencies from the given protocol, falling back to earlier ones."""
        p = str(protocol_number)
        if p in self.data:
            freq = self.data[p].Freq
            if freq is not None and len(freq) > 0:
                return freq
    
        for i in range(int(protocol_number) - 1, -1, -1):
            key = str(i)
            if key in self.data:
                freq = self.data[key].Freq
                if freq is not None and len(freq) > 0:
                    return freq
    
        return np.array([])
    
    def get_last_bvec(self, protocol_number: str) -> Optional[np.ndarray]:
        """Retrieve B_vec from the given protocol, falling back to earlier ones."""
        for i in range(int(protocol_number), -1, -1):
            key = str(i)
            if key in self.data:
                bv = self.data[key].B_vec
                if bv is not None:
                    return bv
        return None


def compute_rotational_constants(conf: 'Conformer', protocol_number: str) -> None:
    """Compute and store principal rotational constants from conformer geometry.

    Calculates the three principal rotational constants (B_a, B_b, B_c) from
    the inertia tensor derived from the conformer's atomic positions and masses.

    Skips if the EnergyRecord for this protocol already has B set, so it is
    safe to call unconditionally after any calculation path (ML, QM parser, or
    checkpoint restore).

    Args:
        conf: Conformer with atoms and last_geometry populated.
        protocol_number: Protocol step number for the EnergyRecord target.

    Returns:
        None

    Raises:
        KeyError: If no EnergyRecord exists for the given protocol number.
    """
    protocol_number = str(protocol_number)
    if protocol_number not in conf.energies:
        raise KeyError(
            f"No EnergyRecord for protocol {protocol_number} in conformer "
            f"{conf.number}. Cannot compute rotational constants."
        )
    record = conf.energies[protocol_number]
    if record.B is not None:
        return

    atoms = Atoms(
        symbols="".join(list(conf.atoms)),
        positions=conf.last_geometry,
    )
    moments = atoms.get_moments_of_inertia()  # [amu·Å²]
    B_vec = np.divide(
        ROT_CONST_FACTOR, moments,
        out=np.zeros_like(moments),
        where=moments > 1e-6,
    )  # [cm⁻¹]

    conf.energies.set(protocol_number, "B", float(np.linalg.norm(B_vec)))
    conf.energies.set(protocol_number, "B_vec", B_vec)


def compute_thermochemistry(
    conf: 'Conformer', protocol_number, energy, freqs,
    temperature, linear, cut_off, alpha, P, mult,
) -> None:
    rec = conf.energies[protocol_number]
    pos_freq = freqs[freqs > 0]
    if len(pos_freq) > 0 and rec.B_vec is not None:
        try:
            g, zpve, h_val, s_val = free_gibbs_energy(
                SCF=energy, T=temperature, freq=pos_freq,
                mw=conf.weight_mass, B=rec.B_vec, m=mult,
                linear=linear, cut_off=cut_off, alpha=alpha, P=P,
            )
            conf.energies.set(protocol_number, "G", g)
            conf.energies.set(protocol_number, "G_E", g - energy)
            conf.energies.set(protocol_number, "zpve", zpve)
            conf.energies.set(protocol_number, "H", h_val)
            conf.energies.set(protocol_number, "S", s_val)
        except Exception:
            pass


def copy_thermochemical_corrections(conf: 'Conformer', target_protocol: str) -> None:
    for p in range(target_protocol - 1, -1, -1):
        if p in conf.energies:
            prev = conf.energies[p]
            if not np.isnan(prev.G_E):
                rec = conf.energies[target_protocol]
                for attr in ("G_E", "zpve", "H", "S"):
                    setattr(rec, attr, getattr(prev, attr))
                rec.G = rec.E + rec.G_E
                return