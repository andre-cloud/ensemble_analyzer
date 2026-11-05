"""Empirical dispersion corrections as a standalone protocol keyword.

The user writes an engine-agnostic name (``"D3BJ"``, ``"D4"``, ...) in the
``dispersion`` field; each QM driver translates it into its own syntax:

* ORCA     : simple-input keyword           ``! D3BJ``
* Gaussian : route keyword                  ``EmpiricalDispersion=GD3BJ``
* NWChem   : directive inside the DFT block ``disp vdw 4``
"""

import re
from typing import Optional


# Canonical name -> engine-specific value (missing key = not supported)
DISPERSION_KEYWORDS: dict[str, dict[str, object]] = {
    "orca": {
        "D2"    : "D2",
        "D3ZERO": "D3ZERO",
        "D3BJ"  : "D3BJ",
        "D4"    : "D4",
    },
    "gaussian": {
        "D2"    : "GD2",
        "D3ZERO": "GD3",
        "D3BJ"  : "GD3BJ",
        "PFD"   : "PFD",
    },
    "nwchem": {
        "D2"    : 2,     # disp vdw 2
        "D3ZERO": 3,     # disp vdw 3
        "D3BJ"  : 4,     # disp vdw 4
    },
}

_ALIASES = {
    "d2"        : "D2",
    "gd2"       : "D2",
    "d3zero"    : "D3ZERO",
    "d3(0)"     : "D3ZERO",
    "gd3"       : "D3ZERO",
    "d3bj"      : "D3BJ",
    "d3(bj)"    : "D3BJ",
    "gd3bj"     : "D3BJ",
    "d4"        : "D4",
    "pfd"       : "PFD",
}

_NONE_VALUES = {"", "none", "no", "off", "false"}

# Functionals that already carry their own dispersion term
_BUILTIN_DISPERSION_FUNCTIONALS = {
    "wb97xd", "b97d", "b97d3", "apfd", "b2plypd", "b2plypd3",
    "wb97x-v", "wb97m-v",
}
_DISP_SUFFIX = re.compile(r"-?d(2|3|3bj|3zero|4)(rev)?$", re.IGNORECASE)

# Patterns that reveal a dispersion already requested through add_input
_ADD_INPUT_PATTERNS = {
    "orca"    : re.compile(r"\b(d2|d3|d3bj|d3zero|d4)\b", re.IGNORECASE),
    "gaussian": re.compile(r"empiricaldispersion", re.IGNORECASE),
    "nwchem"  : re.compile(r"^\s*disp\b", re.IGNORECASE | re.MULTILINE),
}


def normalize_dispersion(value: Optional[str]) -> Optional[str]:
    """Return the canonical dispersion name, or None when disabled.

    Raises:
        ValueError: For unknown names, or the ambiguous ``"D3"`` (ORCA reads
            it as D3BJ, Gaussian's GD3 is zero damping).
    """
    if value is None or value is False:
        return None
    key = str(value).strip().lower()
    if key in _NONE_VALUES:
        return None
    if key == "d3":
        raise ValueError(
            "dispersion='D3' is ambiguous (ORCA: D3BJ, Gaussian GD3: zero damping). "
            "Use 'D3BJ' or 'D3ZERO'."
        )
    if key not in _ALIASES:
        raise ValueError(
            f"dispersion={value!r} not recognised. "
            f"Choose one of {sorted(set(_ALIASES.values()))} or null."
        )
    return _ALIASES[key]


def has_builtin_dispersion(functional: str) -> bool:
    """True if the functional name already includes a dispersion correction."""
    f = functional.strip().lower()
    return f in _BUILTIN_DISPERSION_FUNCTIONALS or bool(_DISP_SUFFIX.search(f))


def dispersion_keyword(dispersion: Optional[str], calculator: str) -> Optional[object]:
    """Engine-specific value for a canonical dispersion name (None if unset)."""
    if not dispersion:
        return None
    return DISPERSION_KEYWORDS[calculator.lower()][dispersion]


def validate_dispersion(dispersion: Optional[str], calculator: str,
                        functional: str, add_input: str = "",
                        builtin: bool = False) -> None:
    """Check that the requested dispersion can be applied without double counting.

    Args:
        builtin: True when the method already includes dispersion for other
            reasons (e.g. composite methods, xTB).

    Raises:
        ValueError: If the calculator does not support the correction, or the
            functional / add_input already provide one.
    """
    if not dispersion:
        return
    calc = calculator.lower()
    supported = DISPERSION_KEYWORDS.get(calc)
    if supported is None:
        raise ValueError(
            f"calculator '{calculator}' does not support the 'dispersion' keyword "
            f"(available for {sorted(DISPERSION_KEYWORDS)})."
        )
    if dispersion not in supported:
        raise ValueError(
            f"dispersion '{dispersion}' is not available in {calculator}. "
            f"Choose one of {sorted(supported)}."
        )
    if builtin or has_builtin_dispersion(functional):
        raise ValueError(
            f"functional '{functional}' already includes a dispersion correction: "
            f"remove 'dispersion' or use the plain functional."
        )
    if add_input and _ADD_INPUT_PATTERNS[calc].search(add_input):
        raise ValueError(
            "dispersion is requested both in 'dispersion' and in 'add_input': keep only one."
        )
