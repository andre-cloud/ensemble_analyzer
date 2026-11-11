"""
ΔSCF (orbital-optimised excited states) for ORCA and NWChem: engine-neutral logic.

A protocol step with ``"deltascf": {...}`` optimises a non-aufbau determinant
on top of a reference ground state of the same multiplicity (S0, T1 or D0).
This module holds everything that does not need a QM engine:

* settings normalisation/validation (``normalize_settings``, ``check_protocols``);
* excitations written with relative labels (``HOMO-n->LUMO+m``) and their
  translation to ORCA ``alphaconf``/``betaconf`` lists and NWChem ``swap`` pairs;
* ``resolve_auto``: excitation taken from the dominant transition of a TD-DFT root;
* ``classify``: which spin purification applies to the target determinant;
* the convergence ladders tried in order for each engine;
* Yamaguchi approximate spin projection (open-shell singlets) and the
  recoupling of mixed determinants for doublets with three open shells
  (Hait et al., J. Chem. Phys. 153, 134108 (2020), Eqs. 1-3).

The validated behaviour (ORCA 6.1.1 / NWChem 7.3.1, formaldehyde and water)
is documented in docs/source/user_guide/deltascf.md.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

EH_TO_EV = 27.211386245981

SUPPORTED_CALCULATORS = ("orca", "nwchem")

DEFAULTS: Dict = {
    "excitation": "auto",      # "auto" | "HOMO->LUMO" | "HOMO-n->LUMO+m"
    "root": 1,                 # TD-DFT root used by "auto" (relative to the step multiplicity)
    "spin": "alpha",           # spin of the excited electron (explicit excitations)
    "alphaconf": None,         # ORCA only: raw alphaconf list (overrides excitation)
    "betaconf": None,          # ORCA only: raw betaconf list (overrides excitation)
    "min_gap_ev": 0.1,         # collapse check vs the reference of the same multiplicity
    "max_attempts": None,      # None = the whole ladder of the engine
    "cation_guess": True,      # last attempts start from the cation orbitals
    "spin_purify": True,       # AP (singlets) / recoupling (doublets, 3 open shells)
    "doublet_branch": "lower", # recoupled doublet stored as E: "lower" | "upper"
    "reference_step": None,    # filled by check_protocols: reused TD-DFT step, or None (#0 computed)
}

_EXCITATION_RE = re.compile(r"^HOMO(?:-(\d+))?->LUMO(?:\+(\d+))?$", re.IGNORECASE)
_FORBIDDEN_IN_ADD_INPUT = ("deltascf", "alphaconf", "betaconf", "max_ovl", "imom", "swap")


# ---------------------------------------------------------------------------
# Excitations
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Excitation:
    """Single excitation HOMO-n -> LUMO+m of one spin, relative to the reference."""
    spin: str                     # "alpha" | "beta"
    n: int = 0                    # HOMO-n
    m: int = 0                    # LUMO+m
    weight: Optional[float] = None
    source: str = "user"          # "user" | "auto" | "raw"

    @property
    def label(self) -> str:
        homo = "HOMO" if self.n == 0 else f"HOMO-{self.n}"
        lumo = "LUMO" if self.m == 0 else f"LUMO+{self.m}"
        return f"{homo}->{lumo} ({self.spin})"


def parse_excitation(text: str) -> Tuple[int, int]:
    """``"HOMO-1->LUMO+2"`` -> (1, 2). Spaces are ignored."""
    match = _EXCITATION_RE.match(re.sub(r"\s+", "", str(text)))
    if not match:
        raise ValueError(
            f"excitation {text!r} not understood: use 'auto', 'HOMO->LUMO' or 'HOMO-n->LUMO+m'"
        )
    return int(match.group(1) or 0), int(match.group(2) or 0)


def orca_conf(n: int, m: int) -> str:
    """ORCA alphaconf/betaconf occupations for HOMO-n -> LUMO+m.

    ORCA lists the occupations starting from HOMO-(Σocc-1) upwards, so the
    excitation is ``[0] + [1]*n + [0]*m + [1]`` (verified against NWChem swaps:
    "0,1" = HOMO->LUMO, "0,1,1" = HOMO-1->LUMO, "0,0,1" = HOMO->LUMO+1).
    """
    return ",".join(str(x) for x in [0] + [1] * n + [0] * m + [1])


def nwchem_swap(n: int, m: int, n_occ: int) -> Tuple[int, int]:
    """Absolute (1-based) NWChem orbital pair for HOMO-n -> LUMO+m, given the
    number of occupied orbitals of that spin in the reference."""
    return n_occ - n, n_occ + 1 + m


def classify(mult: int, exc: Excitation) -> str:
    """Spin purification needed by the target determinant.

    The reference has ``mult - 1`` singly occupied (alpha) orbitals.

    * ``"pure"``: the excitation moves an existing open shell (alpha out of a
      SOMO, or beta into a SOMO): all unpaired electrons stay parallel.
    * ``"singlet_ap"``: closed-shell reference, two antiparallel open shells.
    * ``"doublet_recoupling"``: doublet with three open shells.
    * ``"unsupported"``: four or more open shells (no validated scheme).
    """
    somo = mult - 1
    if (exc.spin == "alpha" and exc.n < somo) or (exc.spin == "beta" and exc.m < somo):
        return "pure"
    open_shells = somo + 2
    return {2: "singlet_ap", 3: "doublet_recoupling"}.get(open_shells, "unsupported")


def expected_s2(mult: int) -> float:
    s = (mult - 1) / 2
    return s * (s + 1)


# ---------------------------------------------------------------------------
# Excitation from a TD-DFT root
# ---------------------------------------------------------------------------

@dataclass
class ExcitedState:
    """One TD-DFT root as parsed from the reference output."""
    root: int
    energy_ev: float
    # (from, to, spin, weight) with engine-native absolute orbital indices
    transitions: List[Tuple[int, int, str, float]]


def resolve_auto(states: Sequence[ExcitedState], root: int, n_alpha: int, n_beta: int,
                 index_base: int) -> Excitation:
    """Dominant transition of ``root`` as a relative excitation.

    Args:
        states: TD-DFT roots of the reference (same multiplicity as the step).
        root: root number (1 = lowest excited state of that multiplicity).
        n_alpha, n_beta: occupied orbitals per spin in the reference.
        index_base: 0 for ORCA (orbitals numbered from 0), 1 for NWChem.
    """
    state = next((s for s in states if s.root == root), None)
    if state is None:
        raise ValueError(
            f"TD-DFT root {root} not found in the reference output "
            f"(roots available: {[s.root for s in states]})"
        )
    if not state.transitions:
        raise ValueError(f"TD-DFT root {root}: no orbital composition printed")
    frm, to, spin, weight = max(state.transitions, key=lambda t: t[3])
    n_occ = n_alpha if spin == "alpha" else n_beta
    homo = n_occ - 1 + index_base
    n, m = homo - frm, to - (homo + 1)
    if n < 0 or m < 0:
        raise ValueError(f"TD-DFT root {root}: dominant transition {frm}->{to} ({spin}) is not occupied->virtual")
    return Excitation(spin=spin, n=n, m=m, weight=weight, source="auto")


# ---------------------------------------------------------------------------
# Spin purification
# ---------------------------------------------------------------------------

def yamaguchi(e_mix: float, s2_mix: float, e_t: float, s2_t: float) -> float:
    """Approximate spin projection of an open-shell singlet (Yamaguchi):
    E_S = E_mix + (E_mix - E_T) * <S²>_mix / (<S²>_T - <S²>_mix).
    With <S²>_mix = 1 and <S²>_T = 2 it reduces to 2 E_mix - E_T."""
    return e_mix + (e_mix - e_t) * s2_mix / (s2_t - s2_mix)


def recouple_doublet(e_q: float, e_m1: float, e_m2: float, e_m3: float) -> Tuple[float, float]:
    """Two spin-pure doublets from the quartet |↑↑↑> and the three mixed
    determinants |↓↑↑>, |↑↓↑>, |↑↑↓> (Hait et al. 2020, Eqs. 2-3).

    Returns:
        (E_lower, E_upper)
    """
    total = e_m1 + e_m2 + e_m3 - e_q
    root = math.sqrt(2 * ((e_m1 - e_m2) ** 2 + (e_m2 - e_m3) ** 2 + (e_m3 - e_m1) ** 2))
    return 0.5 * (total - root), 0.5 * (total + root)


def exchange_couplings(e_q: float, e_m1: float, e_m2: float, e_m3: float) -> Tuple[float, float, float]:
    """K12, K13, K23 from the four determinant energies (must be > 0)."""
    k12 = (e_m1 + e_m2 - e_q - e_m3) / 2
    k13 = (e_m1 + e_m3 - e_q - e_m2) / 2
    k23 = (e_m2 + e_m3 - e_q - e_m1) / 2
    return k12, k13, k23


def doublet_determinants(exc: Excitation) -> Tuple[Dict[str, Tuple[int, List[Excitation]]], str]:
    """Determinants of a doublet excitation with three open shells i (doubly
    occupied in D0), s (SOMO), a (virtual), as single excitations relative to
    the aufbau determinant of the given multiplicity.

    Returns:
        ({name: (mult, [excitations])}, target_name). ``M1`` = β in i (α i→a),
        ``M2`` = β in s, ``M3`` = β in a (β i→a), ``Q`` = quartet |↑↑↑>.
    """
    if exc.spin == "alpha":
        n, m, target = exc.n, exc.m, "M1"
    else:  # β HOMO-n -> LUMO+m (m >= 1, above the SOMO) = α-indexed i = HOMO-(n+1), a = LUMO+(m-1)
        n, m, target = exc.n + 1, exc.m - 1, "M3"
    if n < 1 or m < 0:
        raise ValueError(f"{exc.label} does not open three shells on a doublet")
    q: List[Excitation] = []
    if m >= 1:
        q.append(Excitation("alpha", 0, m - 1))
    if n >= 2:
        q.append(Excitation("beta", n - 2, 0))
    dets = {
        "M1": (2, [Excitation("alpha", n, m)]),
        "M2": (2, [Excitation("alpha", 0, m), Excitation("beta", n - 1, 0)]),
        "M3": (2, [Excitation("beta", n - 1, m + 1)]),
        "Q": (4, q),
    }
    return dets, target


# ---------------------------------------------------------------------------
# Convergence ladders (validated: ORCA 6.1.1, NWChem 7.3.1)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Attempt:
    name: str
    guess: str          # "reference" | "cation"
    keywords: str = ""  # ORCA simple-input keywords
    scf: str = ""       # ORCA %scf options
    dft: Optional[Dict] = None  # NWChem dft block entries


ORCA_LADDER: Tuple[Attempt, ...] = (
    Attempt("GMF", "reference", "GMF", "SOSCFUpdateSPOEst false MaxIter 200"),
    Attempt("GMF, short steps", "reference", "GMF", "SOSCFUpdateSPOEst false SOSCFMaxStep 0.1 MaxIter 300"),
    Attempt("FreezeAndRelease", "reference", "FreezeAndRelease SlowConv", "MaxIter 300"),
    Attempt("FreezeAndRelease, very slow", "reference", "FreezeAndRelease VerySlowConv",
            "PMOM true SOSCFMaxStep 0.05 MaxIter 500"),
    Attempt("GMF from cation", "cation", "GMF", "SOSCFUpdateSPOEst false MaxIter 300"),
    Attempt("FreezeAndRelease from cation", "cation", "FreezeAndRelease VerySlowConv",
            "PMOM true SOSCFMaxStep 0.05 MaxIter 500"),
)

# ORCA 6.0 has DeltaSCF but neither GMF nor FreezeAndRelease
ORCA60_LADDER: Tuple[Attempt, ...] = (
    Attempt("MOM", "reference", "SlowConv", "MaxIter 300"),
    Attempt("PMOM, very slow", "reference", "VerySlowConv", "PMOM true SOSCFMaxStep 0.05 MaxIter 500"),
    Attempt("MOM from cation", "cation", "SlowConv", "MaxIter 300"),
    Attempt("PMOM from cation", "cation", "VerySlowConv", "PMOM true SOSCFMaxStep 0.05 MaxIter 500"),
)

NWCHEM_LADDER: Tuple[Attempt, ...] = (
    Attempt("MOM", "reference", dft={"max_ovl": None}),
    Attempt("IMOM", "reference", dft={"imom": None}),
    Attempt("MOM, damped", "reference", dft={"max_ovl": None, "convergence": {"damp": 70, "ncydp": 50},
                                              "iterations": 300}),
    Attempt("MOM from cation", "cation", dft={"max_ovl": None, "iterations": 300}),
    Attempt("IMOM from cation, damped", "cation", dft={"imom": None, "convergence": {"damp": 70, "ncydp": 50},
                                                        "iterations": 300}),
)


def orca_version() -> Tuple[int, int]:
    """(major, minor) from ORCAVERSION, (0, 0) when unknown."""
    numbers = re.findall(r"\d+", os.getenv("ORCAVERSION", "") or "")
    if not numbers:
        return 0, 0
    return int(numbers[0]), int(numbers[1]) if len(numbers) > 1 else 0


def ladder(calculator: str, settings: Dict) -> Tuple[Attempt, ...]:
    """Attempts tried in order for this engine and these settings."""
    if calculator == "nwchem":
        attempts = NWCHEM_LADDER
    else:
        attempts = ORCA60_LADDER if orca_version() == (6, 0) else ORCA_LADDER
    if not settings.get("cation_guess", True):
        attempts = tuple(a for a in attempts if a.guess != "cation")
    if settings.get("max_attempts"):
        attempts = attempts[: int(settings["max_attempts"])]
    return attempts


def cation_multiplicity(mult: int) -> int:
    """Removing one electron: singlet -> doublet, otherwise one unpaired electron less."""
    return 2 if mult == 1 else mult - 1


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

def normalize_settings(raw, number: int, calculator: str, add_input: str = "",
                       nroots: Optional[int] = None) -> Optional[Dict]:
    """Validate the ``deltascf`` protocol key and fill in the defaults.

    Returns a plain dict (it is serialised with the protocol) or None.
    """
    if raw is None or raw is False:
        return None
    if raw is True:          # "deltascf": true / {} -> all defaults
        raw = {}
    if not isinstance(raw, dict):
        raise ValueError(f"Protocol {number}: 'deltascf' must be a dictionary of options")
    unknown = set(raw) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"Protocol {number}: unknown deltascf option(s) {sorted(unknown)}. "
                         f"Valid: {sorted(k for k in DEFAULTS if k != 'reference_step')}")
    settings = {**DEFAULTS, **raw}

    calc = (calculator or "").lower()
    if calc not in SUPPORTED_CALCULATORS:
        raise ValueError(f"Protocol {number}: deltascf is available for {SUPPORTED_CALCULATORS}, not '{calculator}'")
    if nroots:
        raise ValueError(f"Protocol {number}: a deltascf step cannot also request TD-DFT (nroots); "
                         "put TD-DFT in the reference step")
    lowered = (add_input or "").lower()
    clash = [kw for kw in _FORBIDDEN_IN_ADD_INPUT if re.search(rf"\b{kw}\b", lowered)]
    if clash:
        raise ValueError(f"Protocol {number}: remove {clash} from add_input, deltascf writes them itself")

    exc = str(settings["excitation"]).strip()
    if exc.lower() != "auto":
        parse_excitation(exc)
    settings["excitation"] = "auto" if exc.lower() == "auto" else re.sub(r"\s+", "", exc).upper()
    settings["spin"] = str(settings["spin"]).strip().lower()
    if settings["spin"] not in ("alpha", "beta"):
        raise ValueError(f"Protocol {number}: deltascf spin must be 'alpha' or 'beta'")
    if int(settings["root"]) < 1:
        raise ValueError(f"Protocol {number}: deltascf root must be >= 1")
    settings["root"] = int(settings["root"])
    for key in ("alphaconf", "betaconf"):
        if settings[key] is not None:
            if calc != "orca":
                raise ValueError(f"Protocol {number}: deltascf '{key}' is an ORCA-only override")
            value = str(settings[key]).replace(" ", "")
            if not re.fullmatch(r"[0-2](,[0-2])*", value):
                raise ValueError(f"Protocol {number}: deltascf {key}={settings[key]!r} is not a list of occupations")
            settings[key] = value
    if float(settings["min_gap_ev"]) < 0:
        raise ValueError(f"Protocol {number}: deltascf min_gap_ev must be >= 0")
    settings["min_gap_ev"] = float(settings["min_gap_ev"])
    if settings["max_attempts"] is not None and int(settings["max_attempts"]) < 1:
        raise ValueError(f"Protocol {number}: deltascf max_attempts must be >= 1 (or null)")
    if settings["doublet_branch"] not in ("lower", "upper"):
        raise ValueError(f"Protocol {number}: deltascf doublet_branch must be 'lower' or 'upper'")
    if calc == "orca":
        major, minor = orca_version()
        if 0 < major < 6:
            raise ValueError(f"Protocol {number}: deltascf needs ORCA >= 6.0 (ORCAVERSION={os.getenv('ORCAVERSION')})")
    return settings


_LEVEL_FIELDS = ("calculator", "functional", "basis", "dispersion", "charge", "mult")


def _level(p) -> Tuple:
    solvent = str(p.solvent) if getattr(p, "solvent", None) else None
    return tuple(str(getattr(p, f, "")).lower() for f in _LEVEL_FIELDS) + (solvent,)


def check_protocols(protocols: Sequence, log=None) -> None:
    """Cross-step checks of the deltascf steps, run once all protocols are known.

    ``read_orbitals`` of a deltascf step must point to a TD-DFT step (nroots)
    that is not itself a deltascf step: reading ΔSCF orbitals and applying the
    excitation again collapses to the ground state. A reference at a different
    level of theory is not reused: the reference (#0) is then computed at the
    ΔSCF level. The chosen reference is stored in ``deltascf['reference_step']``.
    """
    by_number = {int(p.number): p for p in protocols}
    for p in protocols:
        settings = getattr(p, "deltascf", None)
        if not settings:
            continue
        ref = str(getattr(p, "read_orbitals", "") or "").strip()
        settings["reference_step"] = None
        if not ref:
            continue
        q = by_number.get(int(ref)) if ref.isdigit() else None
        if q is None:
            raise ValueError(f"Protocol {p.number}: read_orbitals={ref!r} is not a protocol step")
        if getattr(q, "deltascf", None):
            raise ValueError(
                f"Protocol {p.number}: read_orbitals points to the deltascf step {q.number}; "
                "a ΔSCF step must start from ground-state orbitals (TD-DFT reference step)"
            )
        if not getattr(q, "nroots", None):
            raise ValueError(
                f"Protocol {p.number}: the reference step {q.number} (read_orbitals) must be a TD-DFT "
                "step (nroots): it provides ground-state orbitals, E(GS) and the state composition"
            )
        if int(q.number) >= int(p.number):
            raise ValueError(f"Protocol {p.number}: the reference step {q.number} must come before it")
        if _level(q) != _level(p):
            if log is not None:
                log.warning(
                    f"Protocol {p.number}: reference step {q.number} is at a different level of theory; "
                    "the TD-DFT reference (#0) will be computed at the ΔSCF level (one extra TD-DFT)"
                )
            continue
        settings["reference_step"] = int(q.number)
