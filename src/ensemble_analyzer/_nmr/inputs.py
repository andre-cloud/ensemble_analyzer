"""Atom selection and engine keywords for NMR steps.

Protocol indices are 0-based. ORCA ``Nuclei``, Gaussian ``ReadAtoms`` and NWChem ``property``
take 1-based indices (checked on ORCA 6.1, Gaussian 16 A.03, NWChem 7). Syntax notes:

- ORCA: with a list of atoms the isotope must be given (``ist``), otherwise heteronuclear J
  come out as zero; ``all X`` uses the default magnetic isotope, ``ist`` is still written.
- NWChem: spin-spin couplings need an unrestricted wave function (``odft``).
"""

from itertools import combinations
from typing import Dict, List, Sequence

from ensemble_analyzer._nmr.nuclei import element_of, mass_number


def selected_atoms(settings: Dict, symbols: Sequence[str]) -> Dict[str, List[int]]:
    """0-based atoms computed for each requested nucleus (order of ``settings['nuclei']``)."""
    n = len(symbols)
    subset = settings.get("atoms")
    if subset is not None:
        bad = [i for i in subset if i >= n]
        if bad:
            raise ValueError(f"nmr atoms {bad} out of range: the molecule has {n} atoms (0-based)")
    pool = range(n) if subset is None else sorted(set(subset))
    out = {}
    for nuc in settings["nuclei"]:
        el = element_of(nuc)
        out[nuc] = [i for i in pool if str(symbols[i]).capitalize() == el]
    if not any(out.values()):
        raise ValueError(f"nmr: no atom of {settings['nuclei']} among the selected atoms")
    return out


def all_selected(settings: Dict, symbols: Sequence[str]) -> List[int]:
    return sorted(i for idx in selected_atoms(settings, symbols).values() for i in idx)


def orca_block(settings: Dict, symbols: Sequence[str]) -> str:
    """``%eprnmr`` block (written after the coordinates)."""
    flags = []
    if settings["shielding"]:
        flags.append("shift")
    if settings["coupling"]:
        flags.append("ssfc" if settings.get("coupling_terms") == "fc" else "ssall")
    subset = settings.get("atoms") is not None
    lines = ["%eprnmr"]
    for nuc, idx in selected_atoms(settings, symbols).items():
        if not idx:
            continue
        target = ", ".join(str(i + 1) for i in idx) if subset else f"all {element_of(nuc)}"
        lines.append(f"  Nuclei = {target} {{{', '.join(flags)}, ist = {mass_number(nuc)}}}")
    lines.append("end")
    return "\n".join(lines)


def gaussian_keywords(settings: Dict, symbols: Sequence[str]):
    """(route keyword, ReadAtoms line or "")."""
    opts = ["GIAO"]
    atoms_line = ""
    if settings["coupling"]:
        opts.append("SpinSpin")
        if settings.get("coupling_terms") == "fc":
            opts.append("FCOnly")
        chosen = all_selected(settings, symbols)
        if len(chosen) != len(symbols):
            opts.append("ReadAtoms")
            atoms_line = "atoms=" + ",".join(str(i + 1) for i in chosen)
    route = "NMR=GIAO" if len(opts) == 1 else f"NMR=({','.join(opts)})"
    return route, atoms_line


def nwchem_property(settings: Dict, symbols: Sequence[str]) -> Dict:
    """``property`` block: shielding / spinspin, with explicit atoms only for a subset."""
    chosen = all_selected(settings, symbols)
    subset = len(chosen) != len(symbols)
    prop = {}
    if settings["shielding"]:
        prop["shielding"] = f"{len(chosen)} " + " ".join(str(i + 1) for i in chosen) if subset else None
    if settings["coupling"]:
        if subset:
            pairs = list(combinations(chosen, 2))
            prop["spinspin"] = f"{len(pairs)} " + " ".join(f"{a + 1} {b + 1}" for a, b in pairs)
        else:
            prop["spinspin"] = None
    return prop
