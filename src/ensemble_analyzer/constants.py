import os
from pathlib import Path

# --- Light constants (no numpy/scipy needed) ---

DEBUG = bool(os.getenv("DEBUG"))
MAX_TRY = 5

CHIRALS = ["VCD", "ECD"]
GRAPHS = ['IR', 'VCD', 'UV', 'ECD']

CONVERT_B = {'GHz': 29.979000}

VIBRO_OR_ELECTRO = {
    'IR': 'vibro',
    'VCD': 'vibro',
    'UV': 'electro',
    'ECD': 'electro',
}

LOG_FORMAT = "%(message)s"

EH_TO_KCAL = 627.5096080305927
EV_TO_EH = 1.0 / 27.211386245981  # NIST CODATA 2022
CAL_TO_J = 4.186

MARKERS = [
    ".", ",", "o", "v", "^", "<", ">", "1", "2", "3",
    "4", "8", "s", "p", "*", "h", "H", "+", "x", "D",
    "d", "|", "_", "P", "X",
]

MIN_RETENTION_RATE = 0.2
DEFAULT_RESOLUTION = 500
MIN_CONFORMERS_FOR_PCA = 50
MIN_WEIGHTED_VALUE = 0.15

regex_parsing = {
    "orca": {"ext": "out"},
    "gaussian": {"ext": "log"},
    "nwchem": {"ext": "nwo"},
    "tblite": {"ext": None},
    "aimnet": {"ext": None},
    "uma": {"ext": None},
    "mace": {"ext": None},
}

# --- Lazy — initialised on first access to scipy/numpy-dependent constants ---

_LAZY_INIT = False

def _ensure_lazy():
    global _LAZY_INIT, np, R, c, h, electron_volt, Boltzmann, N_A, \
           FACTOR_EV_NM, FACTOR_EV_CM_1, J_TO_H, AMU_TO_KG, ROT_CONST_FACTOR
    if _LAZY_INIT:
        return

    import numpy as _np
    from scipy.constants import R as _R, c as _c, h as _h, \
        electron_volt as _ev, Boltzmann as _B, N_A as _NA
    from scipy.constants import physical_constants

    _J_TO_H = physical_constants["joule-hartree relationship"][0]
    _AMU_TO_KG = physical_constants["atomic mass constant"][0]

    g = globals()
    g['np'] = _np
    g['R'] = _R
    g['h'] = _h
    g['electron_volt'] = _ev
    g['Boltzmann'] = _B
    g['N_A'] = _NA
    g['FACTOR_EV_NM'] = _h * _c / (10 ** -9 * _ev)
    g['FACTOR_EV_CM_1'] = 1 / 8065.544
    g['c'] = _c * 100
    g['J_TO_H'] = _J_TO_H
    g['AMU_TO_KG'] = _AMU_TO_KG
    g['ROT_CONST_FACTOR'] = _h / (8 * np.pi ** 2 * g['c'] * _AMU_TO_KG * 1e-20)

    _LAZY_INIT = True


# --- Functions (their bodies depend on lazy constants) ---


def eV_to_nm(eV: "np.ndarray") -> "np.ndarray":
    _ensure_lazy()
    eV = np.maximum(eV, 1e-2)
    return FACTOR_EV_NM / eV


def boltzmann_distribution(
    energies: "np.ndarray",
    temperature: float,
) -> tuple["np.ndarray", "np.ndarray"]:
    _ensure_lazy()
    if energies.size == 0:
        return np.array([]), np.array([])
    rel_energies = energies - energies.min()
    exponent = -(rel_energies * EH_TO_KCAL * 1000 * CAL_TO_J) / (R * temperature)
    weights = np.exp(exponent)
    population = weights / weights.sum()
    return rel_energies * EH_TO_KCAL, population


def ordinal(n: int) -> str:
    return "%d-%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10:: 4])


def compute_boltzmann_populations(conformers, protocol_number: int, temperature: float) -> None:
    import numpy as np
    active = [c for c in conformers if c.active]
    if not active:
        return
    energies = np.array([c.get_energy(protocol_number=protocol_number) for c in active])
    rel_en, pops = boltzmann_distribution(energies, temperature)
    for c, rel_e, p in zip(active, rel_en, pops):
        c.energies.set(protocol_number, 'Pop', p * 100)
        c.energies.set(protocol_number, 'Erel', rel_e)


def get_models_dir(calculator_name: str, create: bool = True) -> Path:
    env_path = os.getenv("ENAN_MODELS_DIR")
    base_dir = Path(env_path) if env_path else Path.home() / ".ensemble_analyzer" / "models"
    models_dir = base_dir / calculator_name
    if create:
        models_dir.mkdir(parents=True, exist_ok=True)
    return models_dir


def __getattr__(name):
    heavy = {'np', 'R', 'c', 'h', 'electron_volt', 'Boltzmann', 'N_A',
             'FACTOR_EV_NM', 'FACTOR_EV_CM_1', 'J_TO_H', 'AMU_TO_KG',
             'ROT_CONST_FACTOR'}
    if name in heavy:
        _ensure_lazy()
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
