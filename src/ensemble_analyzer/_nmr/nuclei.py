"""Nuclear data for NMR (spin, gyromagnetic ratio, natural abundance, shift range).

Data in ``parameters_file/nuclei.json``, keyed by isotope label (``"1H"``, ``"13C"``, ``"15N"``).
``gamma`` is in rad s^-1 T^-1, ``abundance`` is a fraction, ``range`` a default δ window in ppm.
"""

import json
import re
from importlib.resources import files
from typing import Dict, Tuple

NUCLEUS_DATA: Dict[str, Dict] = json.loads(
    files("ensemble_analyzer").joinpath("parameters_file/nuclei.json").read_text()
)

_MASS_FIRST = re.compile(r"^\s*(\d+)\s*-?\s*([A-Za-z]{1,2})\s*$")
_ELEMENT_FIRST = re.compile(r"^\s*([A-Za-z]{1,2})\s*-?\s*(\d+)\s*$")


def normalize_nucleus_name(label: str) -> str:
    """Canonical isotope label: ``"1h"``, ``"13-C"``, ``"C13"`` → ``"1H"``, ``"13C"``.

    Raises:
        ValueError: not an isotope label, or an isotope missing from ``nuclei.json``.
    """
    text = str(label)
    m = _MASS_FIRST.match(text)
    if m:
        mass, sym = m.groups()
    else:
        m = _ELEMENT_FIRST.match(text)
        if not m:
            raise ValueError(
                f"Invalid nucleus {label!r}: give an isotope such as '1H', '13C', '15N'"
            )
        sym, mass = m.groups()
    key = f"{int(mass)}{sym.capitalize()}"
    if key not in NUCLEUS_DATA:
        raise ValueError(f"Unknown nucleus {label!r}. Known: {', '.join(NUCLEUS_DATA)}")
    return key


def element_of(nucleus: str) -> str:
    return NUCLEUS_DATA[normalize_nucleus_name(nucleus)]["element"]


def mass_number(nucleus: str) -> int:
    return int(re.match(r"\d+", normalize_nucleus_name(nucleus)).group())


def get_nucleus_spin(nucleus: str) -> float:
    return float(NUCLEUS_DATA[normalize_nucleus_name(nucleus)]["spin"])


def gamma(nucleus: str) -> float:
    return float(NUCLEUS_DATA[normalize_nucleus_name(nucleus)]["gamma"])


def abundance(nucleus: str) -> float:
    return float(NUCLEUS_DATA[normalize_nucleus_name(nucleus)]["abundance"])


def shift_range(nucleus: str) -> Tuple[float, float]:
    lo, hi = NUCLEUS_DATA[normalize_nucleus_name(nucleus)]["range"]
    return float(lo), float(hi)


def larmor_mhz(nucleus: str, frequency_1h: float) -> float:
    """Larmor frequency (MHz) of ``nucleus`` on a spectrometer of ``frequency_1h`` MHz for 1H."""
    return float(frequency_1h) * abs(gamma(nucleus) / gamma("1H"))


def shielding_to_shift(sigma: float, sigma_ref: float, delta_ref: float = 0.0) -> float:
    """δ = (σ_ref − σ) / (1 − σ_ref·10⁻⁶) + δ_ref, all in ppm."""
    return (sigma_ref - sigma) / (1.0 - sigma_ref * 1e-6) + delta_ref


def rescale_coupling(j_hz: float, used: Tuple[str, str], wanted: Tuple[str, str]) -> float:
    """J for the ``wanted`` isotope pair from a J computed for the ``used`` pair (J ∝ γA·γB)."""
    return j_hz * gamma(wanted[0]) * gamma(wanted[1]) / (gamma(used[0]) * gamma(used[1]))
