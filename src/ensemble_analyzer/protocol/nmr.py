"""The ``nmr`` protocol key: shieldings/shifts and spin-spin couplings (ORCA, Gaussian, NWChem).

Atom indices in the protocol are 0-based; each engine converts them. See the NMR page of the
documentation for the meaning of every option.
"""

import re
from typing import Dict, List, Optional, Sequence

from ensemble_analyzer._nmr.nuclei import element_of, normalize_nucleus_name
from ensemble_analyzer._nmr.references import (
    RDKIT_HINT,
    bundled_reference,
    needs_rdkit,
    rdkit_available,
)

SUPPORTED_CALCULATORS = ("orca", "gaussian", "nwchem")

DEFAULTS: Dict = {
    "shielding": False,          # isotropic shieldings -> shifts (needs references)
    "coupling": False,           # spin-spin couplings J
    "nuclei": None,              # isotopes, e.g. ["1H", "13C"]
    "atoms": None,               # 0-based subset (None = every atom of the requested nuclei)
    "references": None,          # [{"name", "xyz" | "smiles", "nuclei", "atoms", "delta", "charge", "mult"}]
    "coupling_from": None,       # shielding-only step: NMR step giving J (None = auto, False = none)
    "equivalent": None,          # extra groups of equivalent atoms (0-based)
    "decouple": None,            # {observed: [partners]}
    "coupling_terms": "all",     # "all" | "fc"
    "frequency": 400.0,          # spectrometer 1H frequency, MHz
    "fwhm": 1.0,                 # line width, Hz (number or {nucleus: Hz})
    "opt_from": None,            # step whose level optimises the references (None = auto)
    "reference_opt_step": None,  # filled by check_protocols
    "coupling_step": None,       # filled by check_protocols
}
_INTERNAL_KEYS = ("reference_opt_step", "coupling_step")
REFERENCE_KEYS = ("name", "xyz", "smiles", "nuclei", "atoms", "delta", "charge", "mult")
_ADD_INPUT_CLASH = re.compile(r"(%eprnmr|\bnmr\b|\bspinspin\b|\bproperty\b)", re.IGNORECASE)


def _err(number, msg):
    return ValueError(f"Protocol {number}: nmr {msg}")


def _indices(value, number, what) -> List[int]:
    if not isinstance(value, (list, tuple)) or not all(
        isinstance(i, int) and not isinstance(i, bool) and i >= 0 for i in value
    ):
        raise _err(number, f"'{what}' must be a list of 0-based atom indices")
    return [int(i) for i in value]


def _nuclei(value, number, what="nuclei") -> List[str]:
    if not isinstance(value, (list, tuple)) or not value:
        raise _err(number, f"'{what}' must be a non-empty list of isotopes (e.g. ['1H', '13C'])")
    try:
        out = [normalize_nucleus_name(n) for n in value]
    except ValueError as e:
        raise _err(number, str(e)) from None
    return list(dict.fromkeys(out))


def _reference(raw, number, step_nuclei) -> Dict:
    if not isinstance(raw, dict):
        raise _err(number, "each reference must be a dictionary")
    unknown = set(raw) - set(REFERENCE_KEYS)
    if unknown:
        raise _err(number, f"unknown reference option(s) {sorted(unknown)}. Valid: {list(REFERENCE_KEYS)}")
    if not str(raw.get("name") or "").strip():
        raise _err(number, "every reference needs a 'name'")
    ref = {k: raw.get(k) for k in REFERENCE_KEYS}
    ref["name"] = str(ref["name"]).strip()
    ref["nuclei"] = _nuclei(ref["nuclei"], number, f"references[{ref['name']}].nuclei")
    extra = [n for n in ref["nuclei"] if n not in step_nuclei]
    if extra:
        raise _err(number, f"reference {ref['name']!r}: {extra} not requested in 'nuclei'")
    if not ref["xyz"] and bundled_reference(ref["name"]) is None and not ref["smiles"]:
        raise _err(number, f"reference {ref['name']!r}: give 'xyz' or 'smiles' "
                           "(or the name of a bundled reference such as TMS)")
    if ref["atoms"] is not None:
        ref["atoms"] = _indices(ref["atoms"], number, f"references[{ref['name']}].atoms")
    ref["delta"] = float(ref["delta"] or 0.0)
    return ref


def normalize_settings(raw, number, calculator: str, add_input: str = "", nroots=None,
                       opt: bool = False, freq: bool = False, deltascf=None) -> Optional[Dict]:
    """Validate the ``nmr`` key and fill in the defaults (plain, serialisable dict) or None."""
    if raw is None or raw is False:
        return None
    if not isinstance(raw, dict):
        raise _err(number, "must be a dictionary of options")
    unknown = set(raw) - set(DEFAULTS)
    if unknown:
        raise _err(number, f"unknown option(s) {sorted(unknown)}. "
                           f"Valid: {sorted(k for k in DEFAULTS if k not in _INTERNAL_KEYS)}")
    s = {**DEFAULTS, **raw}

    if (calculator or "").lower() not in SUPPORTED_CALCULATORS:
        raise _err(number, f"is available for {SUPPORTED_CALCULATORS}, not '{calculator}'")
    if opt or freq:
        raise _err(number, "steps must be single points (no opt/freq/ts)")
    if nroots:
        raise _err(number, "steps cannot also run TD-DFT (nroots)")
    if deltascf:
        raise _err(number, "steps cannot be deltascf steps")
    if _ADD_INPUT_CLASH.search(add_input or ""):
        raise _err(number, "keywords are written by EnAn: remove NMR/%eprnmr/property from add_input")

    s["shielding"], s["coupling"] = bool(s["shielding"]), bool(s["coupling"])
    if not (s["shielding"] or s["coupling"]):
        raise _err(number, "needs 'shielding' and/or 'coupling' set to true")
    s["nuclei"] = _nuclei(s["nuclei"], number)
    elements = [element_of(n) for n in s["nuclei"]]
    if len(set(elements)) != len(elements):
        raise _err(number, f"nuclei {s['nuclei']}: only one isotope per element (same element twice)")

    if s["atoms"] is not None:
        s["atoms"] = _indices(s["atoms"], number, "atoms")

    refs = [_reference(r, number, s["nuclei"]) for r in (s["references"] or [])]
    if s["shielding"]:
        if not refs:
            raise _err(number, "shielding needs 'references' (one reference per nucleus)")
        for nuc in s["nuclei"]:
            owners = [r["name"] for r in refs if nuc in r["nuclei"]]
            if not owners:
                raise _err(number, f"{nuc} has no reference")
            if len(owners) > 1:
                raise _err(number, f"{nuc} has more than one reference ({owners})")
    s["references"] = refs

    groups = s["equivalent"] or []
    if not isinstance(groups, list) or any(not isinstance(g, list) or len(g) < 2 for g in groups):
        raise _err(number, "'equivalent' must be a list of groups of at least two atoms")
    s["equivalent"] = [_indices(g, number, "equivalent") for g in groups]

    dec = s["decouple"] or {}
    if not isinstance(dec, dict):
        raise _err(number, "'decouple' must be {observed: [partners]}")
    out = {}
    for obs, partners in dec.items():
        o = _nuclei([obs], number, "decouple")[0]
        if o not in s["nuclei"]:
            raise _err(number, f"decouple: {o} is not in 'nuclei'")
        out[o] = _nuclei(partners, number, "decouple")
    s["decouple"] = out

    s["coupling_terms"] = str(s["coupling_terms"]).strip().lower()
    if s["coupling_terms"] not in ("all", "fc"):
        raise _err(number, "coupling_terms must be 'all' or 'fc'")
    if float(s["frequency"]) <= 0:
        raise _err(number, "frequency (MHz) must be > 0")
    s["frequency"] = float(s["frequency"])
    if isinstance(s["fwhm"], dict):
        fw = {}
        for k, v in s["fwhm"].items():
            if float(v) <= 0:
                raise _err(number, "fwhm (Hz) must be > 0")
            fw[_nuclei([k], number, "fwhm")[0]] = float(v)
        s["fwhm"] = fw
    else:
        if float(s["fwhm"]) <= 0:
            raise _err(number, "fwhm (Hz) must be > 0")
        s["fwhm"] = float(s["fwhm"])

    if s["opt_from"] is not None:
        s["opt_from"] = int(s["opt_from"])
    if s["coupling_from"] is not None and s["coupling_from"] is not False:
        s["coupling_from"] = int(s["coupling_from"])
    return s


def is_ground_state_opt(p) -> bool:
    """Steps whose level can optimise an NMR reference: opt, not TS, ΔSCF or TD-DFT."""
    return bool(getattr(p, "opt", False)) and not getattr(p, "ts", False) \
        and not getattr(p, "deltascf", None) and not getattr(p, "nroots", None)


def _source(ref) -> tuple:
    return (ref.get("xyz"), ref.get("smiles"))


def check_protocols(protocols: Sequence, log=None) -> None:
    """Cross-step resolution, run once all protocols are known.

    Fills ``reference_opt_step`` (level used to optimise the references) and ``coupling_step``
    (step giving J to a shielding-only step); checks reference names and RDKit availability.
    """
    by_number = {int(p.number): p for p in protocols}
    seen_refs: Dict[str, tuple] = {}
    smiles_needed = []
    for p in sorted(protocols, key=lambda q: int(q.number)):
        s = getattr(p, "nmr", None)
        if not s:
            continue
        n = int(p.number)
        earlier = [q for k, q in sorted(by_number.items()) if k < n]

        s["reference_opt_step"] = None
        if s["shielding"]:
            if s.get("opt_from") is not None:
                q = by_number.get(int(s["opt_from"]))
                if q is None or q not in earlier or not is_ground_state_opt(q):
                    raise ValueError(
                        f"Protocol {n}: nmr opt_from={s['opt_from']} must be an earlier optimisation "
                        "step (not TS, ΔSCF or TD-DFT)"
                    )
                s["reference_opt_step"] = int(q.number)
            else:
                cands = [q for q in earlier if is_ground_state_opt(q)]
                if cands:
                    s["reference_opt_step"] = int(cands[-1].number)
                elif log is not None:
                    log.warning(
                        f"Protocol {n}: no optimisation step before this NMR step; the NMR "
                        "references will be optimised at the level of the NMR step"
                    )

        s["coupling_step"] = None
        if s["shielding"] and not s["coupling"] and s.get("coupling_from") is not False:
            def has_j(q):
                qs = getattr(q, "nmr", None)
                return bool(qs and qs["coupling"] and set(qs["nuclei"]) & set(s["nuclei"]))
            if s.get("coupling_from") is not None:
                q = by_number.get(int(s["coupling_from"]))
                if q is None or q not in earlier or not has_j(q):
                    raise ValueError(
                        f"Protocol {n}: nmr coupling_from={s['coupling_from']} must be an earlier NMR "
                        "step with coupling on the same nuclei"
                    )
                s["coupling_step"] = int(q.number)
            else:
                cands = [q for q in earlier if has_j(q)]
                if cands:
                    s["coupling_step"] = int(cands[-1].number)

        for ref in s["references"]:
            src = _source(ref)
            if ref["name"] in seen_refs and seen_refs[ref["name"]] != src:
                raise ValueError(
                    f"Protocol {n}: NMR reference {ref['name']!r} is defined differently in another "
                    "step; use a different name for a different molecule"
                )
            seen_refs[ref["name"]] = src
            if needs_rdkit(ref):
                smiles_needed.append(ref["name"])

    if smiles_needed and not rdkit_available():
        raise ImportError(f"NMR references {sorted(set(smiles_needed))}: {RDKIT_HINT}")
