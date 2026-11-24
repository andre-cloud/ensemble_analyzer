"""Shieldings and spin-spin couplings from ORCA, Gaussian and NWChem outputs.

Every engine is reduced to the isotopes requested in the step:

- ORCA prints each pair with its isotopes and g-factors (``13C gnA = 1.405``); a pair printed
  without an isotope (``1C gnA = 0.000``) has a meaningless J and is dropped.
- Gaussian prints J for its own default isotopes; the reduced coupling K is isotope-independent
  and J = c·K·γA·γB with c = 4.35948e-16 (fitted on Gaussian 16 A.03 output, constant to 1e-5).
- NWChem prints the isotopes (``13-C``); its defaults (e.g. 14N) are rescaled with γ. With an
  atom subset the shielding tensors are numbered 1..n in the order of the subset.
"""

import re
from typing import Dict, Sequence, Tuple

from ensemble_analyzer._nmr.inputs import all_selected, selected_atoms
from ensemble_analyzer._nmr.nuclei import NUCLEUS_DATA, gamma, normalize_nucleus_name
from ensemble_analyzer.conformer.nmr_data import NMRRecord

GAUSSIAN_K_TO_J = 4.35948e-16

_FLOAT = r"-?\d+\.\d+(?:[EeDd][-+]?\d+)?"


def _f(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def _isotope(label: str):
    """Canonical isotope or None for labels such as ``1C`` (no isotope assigned)."""
    try:
        key = normalize_nucleus_name(label)
    except ValueError:
        return None
    return key if key in NUCLEUS_DATA else None


# --- ORCA -------------------------------------------------------------------

def _orca(text: str):
    shield: Dict[int, Tuple[float, float]] = {}
    if "CHEMICAL SHIELDING SUMMARY" in text:
        block = text.split("CHEMICAL SHIELDING SUMMARY")[-1]
        header = block.find("-------  -------")
        rows = block[block.index("\n", header) + 1:].split("\n\n")[0] if header >= 0 else ""
        for m in re.finditer(rf"^\s*(\d+)\s+[A-Za-z]+\s+({_FLOAT})\s+({_FLOAT})\s*$", rows, re.M):
            shield[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    pairs = []
    head = re.compile(
        r"NUCLEUS A =\s*[A-Za-z]+\s*(\d+)\s+NUCLEUS B =\s*[A-Za-z]+\s*(\d+)\s*\n"
        r"\s*\(\s*(\w+)\s+gnA\s*=\s*(-?[\d.]+)\s+(\w+)\s+gnB\s*=\s*(-?[\d.]+)\)"
    )
    for chunk in text.split("NUCLEUS A =")[1:]:
        m = head.match("NUCLEUS A =" + chunk)
        tot = re.search(rf"J\[\d+,\d+\]\(Total\).*?iso=\s*({_FLOAT})", chunk)
        if not (m and tot):
            continue
        a, b, iso_a, gn_a, iso_b, gn_b = m.groups()
        if float(gn_a) == 0.0 or float(gn_b) == 0.0:
            continue
        pairs.append((int(a), int(b), _isotope(iso_a), _isotope(iso_b), float(tot.group(1))))
    return shield, pairs


# --- Gaussian ---------------------------------------------------------------

def _lower_triangle(block: str) -> Dict[Tuple[int, int], float]:
    out, cols = {}, None
    for line in block.splitlines()[1:]:
        parts = line.split()
        if not parts:
            continue
        if all(p.isdigit() for p in parts):
            cols = [int(p) - 1 for p in parts]
            continue
        if not parts[0].isdigit() or cols is None:
            break
        row = int(parts[0]) - 1
        for col, val in zip(cols, parts[1:]):
            if col < row:
                out[(col, row)] = _f(val)
    return out


def _gaussian(text: str):
    shield = {}
    if "Magnetic shielding tensor (ppm)" in text:
        block = text.split("Magnetic shielding tensor (ppm)")[-1].split("spin-spin")[0]
        for m in re.finditer(rf"^\s*(\d+)\s+[A-Za-z]+\s+Isotropic =\s*({_FLOAT})\s+Anisotropy =\s*({_FLOAT})",
                             block, re.M):
            shield[int(m.group(1)) - 1] = (float(m.group(2)), float(m.group(3)))
    reduced = {}
    title = "Total nuclear spin-spin coupling K (Hz):"
    if title in text:
        reduced = _lower_triangle(text.split(title)[-1])
    return shield, reduced


# --- NWChem -----------------------------------------------------------------

def _nwchem(text: str):
    shield = []
    if "Chemical Shielding Tensors" in text:
        block = text.split("Chemical Shielding Tensors")[-1]
        for m in re.finditer(
            rf"Atom:\s+(\d+)\s+([A-Za-z]+).*?isotropic =\s*({_FLOAT}).*?anisotropy =\s*({_FLOAT})",
            block, re.S,
        ):
            shield.append((int(m.group(1)), m.group(2), float(m.group(3)), float(m.group(4))))
    pairs = []
    for m in re.finditer(
        rf"Atom\s+(\d+):\s*(\d+)-([A-Za-z]+)\s+and Atom\s+(\d+):\s*(\d+)-([A-Za-z]+).*?"
        rf"Isotropic Spin-Spin Coupling =\s*({_FLOAT})",
        text, re.S,
    ):
        a, ma, ea, b, mb, eb, j = m.groups()
        pairs.append((int(a) - 1, int(b) - 1, _isotope(ma + ea), _isotope(mb + eb), float(j)))
    return shield, pairs


# --- common -----------------------------------------------------------------

def parse_nmr_output(engine: str, text: str, settings: Dict, symbols: Sequence[str]) -> NMRRecord:
    """NMRRecord (0-based atoms, J for the requested isotopes) of one engine output."""
    engine = engine.lower()
    symbols = [str(s).capitalize() for s in symbols]
    per_nucleus = selected_atoms(settings, symbols)
    wanted = {i: nuc for nuc, idx in per_nucleus.items() for i in idx}

    shield: Dict[int, Tuple[float, float]] = {}
    couplings: Dict[Tuple[int, int], float] = {}

    def add_pair(a, b, iso_a, iso_b, j):
        if a not in wanted or b not in wanted or iso_a is None or iso_b is None:
            return
        j = j * gamma(wanted[a]) * gamma(wanted[b]) / (gamma(iso_a) * gamma(iso_b))
        couplings[(min(a, b), max(a, b))] = j

    if engine == "orca":
        shield, pairs = _orca(text)
        for p in pairs:
            add_pair(*p)
    elif engine == "gaussian":
        shield, reduced = _gaussian(text)
        for (a, b), k in reduced.items():
            if a in wanted and b in wanted and k != 0.0:
                couplings[(a, b)] = GAUSSIAN_K_TO_J * k * gamma(wanted[a]) * gamma(wanted[b])
    elif engine == "nwchem":
        raw, pairs = _nwchem(text)
        chosen = all_selected(settings, symbols)
        renumbered = len(chosen) != len(symbols)
        for n, el, iso, aniso in raw:
            idx = chosen[n - 1] if renumbered and n - 1 < len(chosen) else n - 1
            if idx < len(symbols) and symbols[idx] == el.capitalize():
                shield[idx] = (iso, aniso)
        for p in pairs:
            add_pair(*p)
    else:
        raise ValueError(f"No NMR parser for engine {engine!r}")

    record = NMRRecord()
    if settings["shielding"]:
        record.shielding = {i: v[0] for i, v in shield.items() if i in wanted}
        record.aniso = {i: v[1] for i, v in shield.items() if i in wanted}
    if settings["coupling"]:
        record.coupling = couplings
    return record
