"""Ensemble NMR of one protocol step: averages, equivalent nuclei, shifts, tables and spectra.

Fast exchange: σ of every atom and J of every pair are Boltzmann-averaged over the active
conformers (populations of ``read_population`` or of the step itself), then averaged over the
groups of chemically equivalent nuclei found on the most populated conformer (point group +
topology + rotors, see :mod:`.equivalence`) and over the user ``equivalent`` groups.

Outputs (working directory):

- ``nmr_p<N>.out``: tables (references, atoms, signals, couplings, weights);
- ``nmr_p<N>.json``: the same data for further use;
- ``NMR_<isotope>_p<N>_comp.xy``: spectrum of each nucleus (only with ``shielding``), simulated
  with the exact spin Hamiltonian and the couplings of the step or of ``coupling_step``.
"""

import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from tabulate import tabulate

from ensemble_analyzer._nmr.equivalence import DisjointSet, get_nmr_equivalencies
from ensemble_analyzer._nmr.inputs import selected_atoms
from ensemble_analyzer._nmr.nuclei import element_of, larmor_mhz, shielding_to_shift, shift_range
from ensemble_analyzer._nmr.simulation import SpinGroup, simulate_nucleus, stick_spectrum

J_TABLE_THRESHOLD = 0.1  # Hz


@dataclass
class Group:
    nucleus: str
    members: Tuple[int, ...]
    sigma: Optional[float] = None
    delta: Optional[float] = None


@dataclass
class NMRStepResult:
    protocol_number: int
    symbols: List[str]
    nucleus_of: Dict[int, str]
    sigma: Dict[int, float]
    delta: Dict[int, float]
    j: Dict[Tuple[int, int], float]
    groups: List[Group]
    weights: Dict[int, float]
    references: Dict[str, Dict]
    coupling_step: Optional[int] = None
    spectra: Dict[str, Dict] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Averages
# ---------------------------------------------------------------------------

def _weights(conformers, pop_step) -> Dict[int, float]:
    raw = {}
    for c in conformers:
        rec = c.energies[str(pop_step)] if str(pop_step) in c.energies else None
        pop = getattr(rec, "Pop", np.nan) if rec is not None else np.nan
        raw[c.number] = float(pop) if np.isfinite(pop) else 0.0
    total = sum(raw.values())
    if total <= 0:
        return {n: 1.0 / len(raw) for n in raw} if raw else {}
    return {n: w / total for n, w in raw.items()}


def _groups(conformers, weights, symbols, atoms_with_data, user_groups) -> List[List[int]]:
    best = max(conformers, key=lambda c: weights.get(c.number, 0.0))
    dsu = DisjointSet(len(symbols))
    try:
        for groups in get_nmr_equivalencies(symbols, best.last_geometry).values():
            for g in groups:
                for m in g[1:]:
                    dsu.union(g[0], m)
    except Exception:
        pass  # no automatic equivalence: every atom on its own
    for g in user_groups:
        for m in g[1:]:
            if symbols[m] == symbols[g[0]]:
                dsu.union(g[0], m)
    out = []
    for members in dsu.get_groups().values():
        kept = sorted(i for i in members if i in atoms_with_data)
        if kept:
            out.append(kept)
    return sorted(out)


def average_nmr(conformers: Sequence, protocol, references: Dict[str, Dict]) -> NMRStepResult:
    settings = protocol.nmr
    number = int(protocol.number)
    active = [c for c in conformers if c.active]
    pop_step = protocol.read_population or number
    weights = _weights(active, pop_step)
    symbols = [str(s).capitalize() for s in active[0].atoms] if active else []

    nucleus_of = {i: nuc for nuc, idx in selected_atoms(settings, symbols).items() for i in idx}

    sigma: Dict[int, float] = {}
    if settings["shielding"]:
        acc: Dict[int, List[float]] = {}
        for c in active:
            rec = c.nmr_data.get(number)
            if rec is None:
                continue
            for i, s in rec.shielding.items():
                acc.setdefault(i, [0.0, 0.0])
                acc[i][0] += weights[c.number] * s
                acc[i][1] += weights[c.number]
        sigma = {i: v[0] / v[1] for i, v in acc.items() if v[1] > 0}

    j_step = number if settings["coupling"] else settings.get("coupling_step")
    j: Dict[Tuple[int, int], float] = {}
    if j_step is not None:
        with_j = [c for c in active if c.nmr_data.get(j_step) is not None]
        total = sum(weights[c.number] for c in with_j)
        for c in with_j:
            for pair, val in c.nmr_data.get(j_step).coupling.items():
                j[pair] = j.get(pair, 0.0) + weights[c.number] * val / total

    j_nuclei = settings.get("coupling_nuclei") or (settings["nuclei"] if settings["coupling"] else [])
    by_element = {element_of(n): n for n in list(j_nuclei) + list(settings["nuclei"])}
    for a, b in j:
        for i in (a, b):
            if i not in nucleus_of and symbols[i] in by_element:
                nucleus_of[i] = by_element[symbols[i]]

    delta = {}
    for i, s in sigma.items():
        ref = references.get(nucleus_of.get(i))
        if ref:
            delta[i] = shielding_to_shift(s, ref["sigma"], ref.get("delta", 0.0))

    data_atoms = set(nucleus_of) & (set(sigma) | {i for p in j for i in p} | set(nucleus_of))
    groups = []
    for members in _groups(active, weights, symbols, data_atoms, settings.get("equivalent") or []):
        sig = [sigma[i] for i in members if i in sigma]
        dl = [delta[i] for i in members if i in delta]
        groups.append(Group(nucleus_of[members[0]], tuple(members),
                            float(np.mean(sig)) if sig else None, float(np.mean(dl)) if dl else None))

    return NMRStepResult(number, symbols, nucleus_of, sigma, delta, j, groups, weights,
                         dict(references), j_step)


# ---------------------------------------------------------------------------
# Spectra
# ---------------------------------------------------------------------------

def _group_j(result: NMRStepResult) -> np.ndarray:
    n = len(result.groups)
    J = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            vals = [result.j.get((min(x, y), max(x, y)), 0.0)
                    for x in result.groups[a].members for y in result.groups[b].members]
            J[a, b] = J[b, a] = float(np.mean(vals))
    return J


def _lorentzian_on_grid(x: np.ndarray, centers: np.ndarray, heights: np.ndarray, fwhm: float) -> np.ndarray:
    y = np.zeros_like(x)
    step = x[1] - x[0] if len(x) > 1 else 1.0
    half = int(60 * fwhm / step) + 1
    for c, h in zip(centers, heights):
        k = int(round((c - x[0]) / step))
        lo, hi = max(k - half, 0), min(k + half + 1, len(x))
        if lo >= hi:
            continue
        d = x[lo:hi] - c
        y[lo:hi] += h * fwhm ** 2 / (fwhm ** 2 + 4 * d ** 2)
    return y


def _similarity(xc, yc, xe, ye) -> float:
    yi = np.interp(xe, xc, yc, left=0.0, right=0.0)
    den = np.linalg.norm(yi) * np.linalg.norm(ye)
    return float(100.0 * np.dot(yi, ye) / den) if den > 0 else 0.0


def _fit_shift(x, y, xe, ye, window: float) -> Tuple[float, float]:
    """Additive shift (ppm) maximising the overlap with the experiment (FFT cross-correlation
    on the uniform grid of the computed spectrum), limited to ±window."""
    from scipy.signal import correlate

    step = x[1] - x[0]
    ye_grid = np.interp(x, xe, ye, left=0.0, right=0.0)
    corr = correlate(ye_grid, y, mode="full", method="fft")
    lags = (np.arange(corr.size) - (len(y) - 1)) * step
    corr[np.abs(lags) > window] = -np.inf
    shift = float(lags[int(np.argmax(corr))])
    return shift, _similarity(x + shift, y, xe, ye)


def _spectrum(result: NMRStepResult, settings: Dict, nucleus: str, log) -> Optional[Dict]:
    groups = [SpinGroup(g.nucleus, g.delta if g.delta is not None else 0.0, len(g.members), g.members)
              for g in result.groups]
    if not any(g.nucleus == nucleus and result.groups[k].delta is not None for k, g in enumerate(groups)):
        return None
    freq = settings["frequency"]
    fwhm_hz = settings["fwhm"].get(nucleus, 1.0) if isinstance(settings["fwhm"], dict) else settings["fwhm"]
    nu0 = larmor_mhz(nucleus, freq)
    tr = simulate_nucleus(nucleus, groups, _group_j(result), freq,
                          decouple=settings.get("decouple", {}).get(nucleus, []))
    centers, heights = stick_spectrum(tr, nucleus, freq)
    fwhm = fwhm_hz / nu0

    lo_r, hi_r = shift_range(nucleus)
    pad = max(0.5, 0.01 * (hi_r - lo_r))
    lo, hi = centers.min() - pad, centers.max() + pad
    exp_file = f"nmr_{nucleus.lower()}_ref.dat"
    exp = None
    try:
        exp = np.loadtxt(exp_file)
        lo, hi = min(lo, exp[:, 0].min()), max(hi, exp[:, 0].max())
    except OSError:
        exp = None
    step = fwhm / 5
    n = int(np.ceil((hi - lo) / step)) + 1
    if n > 2_000_000:
        step, n = (hi - lo) / 2_000_000, 2_000_001
    x = lo + step * np.arange(n)
    y = _lorentzian_on_grid(x, centers, heights, fwhm)
    y /= y.max()

    info = {"file": f"NMR_{nucleus}_p{result.protocol_number}_comp.xy", "fwhm_hz": fwhm_hz,
            "frequency_mhz": nu0, "lines": int(len(centers)), "shift_ppm": 0.0, "similarity": None}
    if exp is not None:
        xe, ye = exp[:, 0], exp[:, 1] / np.max(np.abs(exp[:, 1]))
        shift, sim = _fit_shift(x, y, xe, ye, window=0.02 * (hi_r - lo_r))
        x = x + shift
        info.update(shift_ppm=shift, similarity=sim, experimental=exp_file)
        log.info(f"NMR {nucleus} protocol {result.protocol_number}: shift {shift:+.3f} ppm, "
                 f"similarity {sim:.1f}% vs {exp_file}")
    np.savetxt(info["file"], np.column_stack((x, y)), delimiter=" ")
    return info


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def _label(result, i):
    return f"{result.symbols[i]}{i}"


def write_tables(result: NMRStepResult, protocol) -> str:
    n = result.protocol_number
    lines = [f"NMR — protocol {n}: {protocol.functional}/{protocol.basis}"
             + (f" [{protocol.solvent}]" if protocol.solvent else ""),
             "Atom indices are 0-based. Shieldings σ and shifts δ in ppm, couplings J in Hz.", ""]
    if result.references:
        rows = [[nuc, r["reference"], f"{r['sigma']:.3f}", f"{r.get('delta', 0.0):.3f}", r["n_atoms"],
                 r.get("opt_step")] for nuc, r in result.references.items()]
        lines += ["References", tabulate(rows, ["nucleus", "reference", "σ_ref", "δ_ref", "atoms",
                                                "opt. level of protocol"], tablefmt="simple"), ""]
    group_of = {i: k for k, g in enumerate(result.groups) for i in g.members}
    if result.sigma:
        rows = [[i, _label(result, i), result.nucleus_of.get(i, ""), f"{s:.3f}",
                 f"{result.delta[i]:.3f}" if i in result.delta else "", group_of.get(i, "")]
                for i, s in sorted(result.sigma.items())]
        lines += ["Atoms", tabulate(rows, ["index", "atom", "nucleus", "σ", "δ", "group"], tablefmt="simple"), ""]
        rows = [[k, g.nucleus, len(g.members), " ".join(map(str, g.members)),
                 f"{g.sigma:.3f}" if g.sigma is not None else "", f"{g.delta:.3f}" if g.delta is not None else ""]
                for k, g in enumerate(result.groups) if g.sigma is not None]
        lines += ["Signals (equivalent nuclei averaged)",
                  tabulate(rows, ["group", "nucleus", "n", "atoms", "σ", "δ"], tablefmt="simple"), ""]
    if result.j:
        rows = [[_label(result, a), _label(result, b), f"{v:.2f}"]
                for (a, b), v in sorted(result.j.items()) if abs(v) >= J_TABLE_THRESHOLD]
        src = "" if result.coupling_step in (None, n) else f" (from protocol {result.coupling_step})"
        lines += [f"Couplings J{src}, |J| ≥ {J_TABLE_THRESHOLD} Hz",
                  tabulate(rows, ["atom A", "atom B", "J"], tablefmt="simple"), ""]
    if result.spectra:
        rows = [[nuc, s["file"], f"{s['frequency_mhz']:.2f}", s["fwhm_hz"], f"{s['shift_ppm']:+.3f}",
                 "" if s["similarity"] is None else f"{s['similarity']:.1f}"] for nuc, s in result.spectra.items()]
        lines += ["Spectra", tabulate(rows, ["nucleus", "file", "MHz", "fwhm (Hz)", "shift (ppm)", "similarity %"],
                                      tablefmt="simple"), ""]
    rows = [[c, f"{100 * w:.2f}"] for c, w in sorted(result.weights.items())]
    lines += ["Conformer weights", tabulate(rows, ["conformer", "weight %"], tablefmt="simple"), ""]
    text = "\n".join(lines)
    with open(f"nmr_p{n}.out", "w", encoding="utf-8") as f:
        f.write(text)

    data = {
        "protocol": n,
        "atoms": [{"index": i, "element": result.symbols[i], "nucleus": result.nucleus_of.get(i),
                   "sigma": result.sigma.get(i), "delta": result.delta.get(i), "group": group_of.get(i)}
                  for i in sorted(set(result.sigma) | set(result.nucleus_of))],
        "groups": [asdict(g) for g in result.groups],
        "couplings": [[a, b, v] for (a, b), v in sorted(result.j.items())],
        "coupling_step": result.coupling_step,
        "references": result.references,
        "weights": {str(k): v for k, v in result.weights.items()},
        "spectra": result.spectra,
    }
    with open(f"nmr_p{n}.json", "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return f"nmr_p{n}.out"


# ---------------------------------------------------------------------------

def process_nmr_step(conformers: Sequence, protocol, references: Dict[str, Dict], log) -> NMRStepResult:
    """Average, write the tables and (with shielding) the spectra of an NMR step."""
    result = average_nmr(conformers, protocol, references)
    settings = protocol.nmr
    if settings["shielding"]:
        for nuc in settings["nuclei"]:
            try:
                info = _spectrum(result, settings, nuc, log)
            except Exception as e:  # a spectrum must never stop the workflow
                log.warning(f"NMR {nuc} spectrum of protocol {protocol.number} failed: {e}")
                info = None
            if info:
                result.spectra[nuc] = info
    fname = write_tables(result, protocol)
    log.info(f"NMR results of protocol {protocol.number} written to {fname}")
    return result
