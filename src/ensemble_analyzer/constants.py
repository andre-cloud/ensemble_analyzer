from scipy.constants import R, c, h, electron_volt, Boltzmann, N_A
from scipy.constants import physical_constants
import numpy as np
from pathlib import Path

import os

DEBUG = bool(os.getenv("DEBUG"))
MAX_TRY = 5


# Physical CONSTANTS

# R = 8.314462618 J/(mol K)
# h = 6.62607015e-34 J*s
# c = 2.9979245800E+10 cm/s
# Boltzmann = 1.380649e-23 J/K
# J_TO_H = 2.2937122783963e+17 Eh/J
# AMU_TO_KG = 1.6605390666e-27 kg*mol/g

FACTOR_EV_NM = h * c / (10**-9 * electron_volt)
FACTOR_EV_CM_1 = 1 / 8065.544  # to yield eV


def eV_to_nm(eV: np.ndarray) -> np.ndarray:
    """Convert energy in eV to wavelength in nm."""
    eV = np.maximum(eV, 1e-2)
    return FACTOR_EV_NM / eV


c = c * 100  # convert speed of light in cm/s
J_TO_H = physical_constants["joule-hartree relationship"][0]
AMU_TO_KG = physical_constants["atomic mass constant"][0]

EH_TO_KCAL = 627.5096080305927
CAL_TO_J = 4.186


CHIRALS = ["VCD", "ECD"]

GRAPHS = ['IR', 'VCD', 'UV', 'ECD']


def boltzmann_distribution(
    energies: np.ndarray,
    temperature: float,
) -> tuple[np.ndarray, np.ndarray]:
    rel_energies = energies - energies.min()
    exponent = -(rel_energies * EH_TO_KCAL * 1000 * CAL_TO_J) / (R * temperature)
    weights = np.exp(exponent)
    population = weights / weights.sum()
    return rel_energies * EH_TO_KCAL, population

CONVERT_B = {
    'GHz': 29.979000,
}

ROT_CONST_FACTOR = h / (8 * np.pi**2 * c * AMU_TO_KG * 1e-20)

VIBRO_OR_ELECTRO = {
    'IR': 'vibro',
    'VCD': 'vibro',
    'UV': 'electro',
    'ECD': 'electro',
}


# Logger constants
LOG_FORMAT = "%(message)s"


def ordinal(n: int) -> str:
    """Return the ordinal suffix string for an integer."""
    return "%d-%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10:: 4])


regex_parsing = {
    "orca": {"ext": "out"},
    "gaussian": {"ext": "log"},
    "nwchem": {"ext": "nwo"},
    "tblite": {"ext": None},
    "aimnet": {"ext": None},
}

MARKERS = [
    ".", ",", "o", "v", "^", "<", ">", "1", "2", "3",
    "4", "8", "s", "p", "*", "h", "H", "+", "x", "D",
    "d", "|", "_", "P", "X",
]

MIN_RETENTION_RATE = 0.2        # Minimum retention rate
DEFAULT_RESOLUTION = 500        # Grid resolution for contour plots
MIN_CONFORMERS_FOR_PCA = 50     # Minimum conformers needed for meaningful PCA
MIN_WEIGHTED_VALUE = 0.15       # Minimum value for the weight of the autoconvolution


def compute_boltzmann_populations(conformers, protocol_number: int, temperature: float) -> None:
    active = [c for c in conformers if c.active]
    if not active:
        return
    energies = np.array([c.get_energy(protocol_number=protocol_number) for c in active])
    rel_en, pops = boltzmann_distribution(energies, temperature)
    for c, rel_e, p in zip(active, rel_en, pops):
        c.energies.set(protocol_number, 'Pop', p * 100)
        c.energies.set(protocol_number, 'Erel', rel_e)


def get_models_dir(calculator_name: str, create: bool = True) -> Path:
    """
    Return the path to ML model weights for a given calculator.

    Priority:
    1. ``ENAN_MODELS_DIR`` environment variable (shared/cluster setup).
    2. ``~/.ensemble_analyzer/models/`` (local default).

    Args:
        calculator_name (str): Subdirectory name (e.g. ``"aimnet"``).
        create (bool): Create the directory if it does not exist.  Defaults to True.

    Returns:
        Path: Absolute path to the calculator-specific model directory.
    """
    env_path = os.getenv("ENAN_MODELS_DIR")
    base_dir = Path(env_path) if env_path else Path.home() / ".ensemble_analyzer" / "models"
    models_dir = base_dir / calculator_name
    if create:
        models_dir.mkdir(parents=True, exist_ok=True)
    return models_dir