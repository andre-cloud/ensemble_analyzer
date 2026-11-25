"""
ΔSCF runner: reference, convergence ladder, acceptance checks and spin
purification for one conformer of a protocol step with ``deltascf``.

Every sub-calculation starts from the ASE calculator that ``OrcaCalc`` /
``NWChemCalc`` build for the step (same functional, basis, dispersion, solvent,
parallel settings and ``add_input``); only the ΔSCF keywords, charge,
multiplicity, initial orbitals and working directory change. Sub-calculations
live in ``<conf>/protocol_<n>/deltascf/``; the accepted ΔSCF output is copied
to the usual ``<conf>_p<n>_<label>.<ext>`` and parsed by get_conf_parameters
(with the spin-purified energy when one is computed).

Engine-neutral logic and the validation behind every choice are in
``ensemble_analyzer.protocol.deltascf``.
"""

from __future__ import annotations

import copy
import json
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from ase import Atoms

from ensemble_analyzer.protocol import deltascf as ds
from ensemble_analyzer.protocol.deltascf import Attempt, Excitation
from ensemble_analyzer._logger.logger import Logger

# Parser log of the trial runs: their warnings are reported by the runner itself
_QUIET = Logger("enan.deltascf.parser")
_QUIET.disabled = True


@dataclass
class _Run:
    """Outcome of one QM sub-calculation."""
    ok: bool
    output: Optional[Path] = None
    orbitals: Optional[Path] = None
    parser: object = None
    energy: Optional[float] = None
    s2: Optional[float] = None
    reason: Optional[str] = None


@dataclass
class _Reference:
    output: Path
    orbitals: Path
    energy: float
    n_alpha: int
    n_beta: int
    states: list = field(default_factory=list)
    reused_step: Optional[int] = None


@dataclass
class _Target:
    """Occupation change of a determinant w.r.t. the aufbau one of its multiplicity."""
    excitations: List[Excitation] = field(default_factory=list)
    alphaconf: Optional[str] = None   # ORCA raw override
    betaconf: Optional[str] = None


# ---------------------------------------------------------------------------
# Engine adapters
# ---------------------------------------------------------------------------

class _Engine:
    """Builds and runs sub-calculations from the step calculator."""

    name = ""
    index_base = 0
    ext = ""

    def __init__(self, protocol, conf, cpu: int):
        self.protocol = protocol
        self.conf = conf
        self.cpu = cpu

    def parser(self, output: Path):
        """Parser of an engine output; its warnings are reported by the runner."""
        from ensemble_analyzer.parsers import PARSER_REGISTRY
        return PARSER_REGISTRY[self.name](str(output), _QUIET)

    def _variant(self, *, mode: str, charge: int, mult: int, nroots: Optional[int]):
        p = copy.copy(self.protocol)          # no __post_init__: keeps Solvent etc.
        p.read_orbitals = ""
        p.deltascf = None
        p.nmr = None
        p.nroots = nroots
        p.tda = True if nroots else self.protocol.tda
        p.charge, p.mult = charge, mult
        p.opt = bool(self.protocol.opt) if mode == "protocol" else False
        p.freq = bool(self.protocol.freq) if mode == "protocol" else False
        p.ts = bool(self.protocol.ts) if mode == "protocol" else False
        return p

    def _base_calc(self, p):
        return p.get_calculator(cpu=self.cpu, conf=self.conf)

    def build(self, run_dir: Path, *, mode: str = "sp", charge: Optional[int] = None,
              mult: Optional[int] = None, guess: Optional[Path] = None,
              target: Optional[_Target] = None, attempt: Optional[Attempt] = None,
              n_occ: Tuple[int, int] = (0, 0), nroots: Optional[int] = None):
        raise NotImplementedError

    def run(self, calc, geometry: np.ndarray, opt: bool) -> _Run:
        raise NotImplementedError


class _OrcaEngine(_Engine):
    name, index_base, ext, orbital_file = "orca", 0, "out", "orca.gbw"

    def build(self, run_dir, *, mode="sp", charge=None, mult=None, guess=None, target=None,
              attempt=None, n_occ=(0, 0), nroots=None):
        charge = self.protocol.charge if charge is None else charge
        mult = self.protocol.mult if mult is None else mult
        calc, _ = self._base_calc(self._variant(mode=mode, charge=charge, mult=mult, nroots=nroots))
        calc.directory = Path(run_dir)
        simple, blocks = [], []
        confs = self._confs(target) if target is not None else ""
        if mult > 1 or target is not None:
            simple.append("UKS")
        if confs:
            # a target without occupation changes (e.g. the recoupling quartet) is a
            # plain UKS on the reference orbitals: ORCA's aufbau keeps that state
            simple += ["DeltaSCF", attempt.keywords if attempt else ""]
            scf = attempt.scf if attempt else ""
            blocks.append(f"%scf {confs} ConvForced true {scf} end")
        if guess is not None:
            simple.append("MORead")
            blocks.insert(0, f'%moinp "{Path(guess).resolve()}"')
        if simple:
            calc.parameters["orcasimpleinput"] += "\n! " + " ".join(s for s in simple if s) + "\n"
        if blocks:
            calc.parameters["orcablocks"] += "\n" + "\n".join(blocks) + "\n"
        return calc

    @staticmethod
    def _confs(target: _Target) -> str:
        alpha, beta = target.alphaconf, target.betaconf
        for exc in target.excitations:
            conf = ds.orca_conf(exc.n, exc.m)
            if exc.spin == "alpha":
                if alpha:
                    raise ValueError("ORCA allows one alpha excitation per determinant")
                alpha = conf
            else:
                if beta:
                    raise ValueError("ORCA allows one beta excitation per determinant")
                beta = conf
        parts = []
        if alpha:
            parts.append(f"alphaconf {alpha}")
        if beta:
            parts.append(f"betaconf {beta}")
        return " ".join(parts)

    def run(self, calc, geometry, opt):
        run_dir = Path(calc.directory)
        run_dir.mkdir(parents=True, exist_ok=True)
        atoms = Atoms(symbols="".join(self.conf.atoms), positions=geometry)
        try:
            calc.write_inputfiles(atoms, ["energy"])
            calc.template.execute(calc.directory, calc.profile)
        except Exception as exc:  # ORCA usually exits 0 even on errors: the text decides
            error = str(exc)
        else:
            error = None
        output = run_dir / calc.template.outputname
        if not output.exists():
            return _Run(False, reason=error or "no output written")
        return _parse(self, output, run_dir / "orca.gbw", opt)


class _NWChemEngine(_Engine):
    name, index_base, ext, orbital_file = "nwchem", 1, "nwo", None

    def build(self, run_dir, *, mode="sp", charge=None, mult=None, guess=None, target=None,
              attempt=None, n_occ=(0, 0), nroots=None):
        charge = self.protocol.charge if charge is None else charge
        mult = self.protocol.mult if mult is None else mult
        calc, _ = self._base_calc(self._variant(mode=mode, charge=charge, mult=mult, nroots=nroots))
        calc.directory = str(Path(run_dir))
        if charge:
            calc.parameters["charge"] = charge
        else:
            calc.parameters.pop("charge", None)
        dft = calc.parameters.setdefault("dft", {})
        dft["mult"] = mult
        out_vectors = Path(run_dir).resolve() / "deltascf.movecs"
        vectors = f'input "{Path(guess).resolve()}"' if guess is not None else ""
        if target is not None:
            dft["odft"] = None
            swaps = {"alpha": [], "beta": []}
            for exc in target.excitations:
                occ = n_occ[0] if exc.spin == "alpha" else n_occ[1]
                swaps[exc.spin] += list(ds.nwchem_swap(exc.n, exc.m, occ))
            for spin, pairs in swaps.items():
                if pairs:
                    vectors += f" swap {spin} " + " ".join(str(i) for i in pairs)
            entries = dict(attempt.dft or {}) if attempt else {}
            convergence = {"lshift": 0.0, **entries.pop("convergence", {})}
            dft["convergence"] = convergence
            dft.update(entries)
        elif mult > 1:
            dft["odft"] = None
        dft["vectors"] = f'{vectors} output "{out_vectors}"'.strip()
        return calc

    def run(self, calc, geometry, opt):
        run_dir = Path(calc.directory)
        run_dir.mkdir(parents=True, exist_ok=True)
        atoms = Atoms(symbols="".join(self.conf.atoms), positions=geometry)
        error = None
        try:
            calc.write_input(atoms, properties=["energy"])
            calc.execute()
        except Exception as exc:  # non-zero exit code: NWChem failed
            error = str(exc) or "NWChem failed"
        output = run_dir / f"{calc.prefix}.nwo"
        if not output.exists():
            return _Run(False, reason=error or "no output written")
        result = _parse(self, output, run_dir / "deltascf.movecs", opt)
        if error and result.ok:
            result.ok, result.reason = False, "NWChem exited with an error"
        return result


def _parse(engine: _Engine, output: Path, orbitals: Path, opt: bool) -> _Run:
    parser = engine.parser(output)
    reason = parser.deltascf_status(opt=opt)
    energy = parser.parse_energy() or None
    if reason is None and energy is None:
        reason = "no energy in the output"
    return _Run(reason is None, output, orbitals if orbitals.exists() else None, parser,
                energy, parser.parse_s2(), reason)


ENGINES = {"orca": _OrcaEngine, "nwchem": _NWChemEngine}


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

class DeltaSCFRunner:
    """Runs a ΔSCF step for one conformer (called by CalculationExecutor)."""

    def __init__(self, executor):
        self.executor = executor
        self.log = executor.logger

    # -- entry point -------------------------------------------------------

    def run(self, conf, protocol, cpu: int) -> bool:
        start = time.perf_counter()
        settings = protocol.deltascf
        engine = ENGINES[protocol.calculator.lower()](protocol, conf, cpu)
        base = Path(conf.folder) / f"protocol_{protocol.number}" / "deltascf"
        tag = f"Conf {conf.number} ΔSCF"
        geometry = np.asarray(conf.last_geometry, dtype=float)

        inherited = None
        try:
            if settings["excitation"] == "inherit":
                inherited = self._load_inherited(conf, settings["inherit_step"])
            ref = self._reference(engine, settings, conf, protocol, base, geometry)
            target, exc, kind = self._target(settings, ref, protocol.mult, engine, inherited)
        except Exception as err:
            self.log.error(f"{tag}: {err}")
            conf.active = False
            return False

        self.log.info(
            f"{tag}: reference {'step ' + str(ref.reused_step) if ref.reused_step is not None else '#0 (computed)'} "
            f"E = {ref.energy:.6f} Eh; excitation {exc.label if exc else 'raw alphaconf/betaconf'}"
            + (f" (weight {exc.weight:.2f})" if exc and exc.weight is not None else "")
        )
        if exc is not None and exc.weight is not None and exc.weight < 0.6:
            self.log.warning(
                f"{tag}: dominant transition weight {exc.weight:.2f} < 0.6: the state is multiconfigurational, "
                "a single-determinant ΔSCF describes it poorly"
            )

        accepted = self._ladder(engine, settings, protocol, ref, target, base, geometry, kind, tag)
        if accepted is None:
            conf.active = False
            self.log.warning(f"{tag}: no attempt converged to an accepted excited state – deactivating")
            return False

        gap = (accepted.energy - ref.energy) * ds.EH_TO_EV
        if inherited is not None:
            self._compare_inherited_gap(gap, inherited, settings, tag)
        self._save_excitation(base, exc, gap, protocol, accepted.energy)

        energy = accepted.energy
        final_geometry = geometry
        if protocol.opt:
            parsed = accepted.parser.parse_geom()
            if parsed.size:
                final_geometry = parsed
        if settings.get("spin_purify", True):
            energy = self._purify(engine, settings, protocol, ref, exc, kind, accepted,
                                  base, final_geometry, tag)
        elif kind == "singlet_ap":
            self.log.warning(
                f"{tag}: spin purification disabled: the open-shell singlet energy is contaminated "
                "by the triplet (error ~ ΔE(S-T)/2)"
            )

        return self._finalize(conf, protocol, engine, accepted, energy, time.perf_counter() - start)

    # -- reference (#0) ----------------------------------------------------

    def _reference(self, engine: _Engine, settings, conf, protocol, base: Path, geometry) -> _Reference:
        step = settings.get("reference_step")
        if step is not None:
            folder = Path(conf.folder) / f"protocol_{step}"
            if engine.name == "orca":
                output, orbitals = folder / f"{conf.number}_p{step}_orca.out", folder / "orca.gbw"
            else:
                label = f"{conf.number}_p{step}_nwchem"
                output, orbitals = folder / f"{label}.nwo", folder / label / f"{label}.movecs"
            if output.exists() and orbitals.exists():
                return self._reference_from(engine, output, orbitals, protocol.mult, step)
            self.log.warning(
                f"Conf {conf.number}: output/orbitals of the reference step {step} not found, computing #0"
            )
        # "inherit": the excitation is known, a ground-state SCF is enough
        # (guess orbitals, E(GS), occupied orbitals for NWChem swaps)
        nroots = None if settings["excitation"] == "inherit" else max(int(settings.get("root", 1)) + 2, 5)
        calc = engine.build(base / "reference", mode="sp", nroots=nroots)
        result = engine.run(calc, geometry, opt=False)
        if not result.ok or result.orbitals is None:
            raise RuntimeError(f"reference GS+TD-DFT (#0) failed: {result.reason}")
        return self._reference_from(engine, result.output, result.orbitals, protocol.mult, None)

    @staticmethod
    def _reference_from(engine, output: Path, orbitals: Path, mult: int, step) -> _Reference:
        parser = engine.parser(output)
        electrons = parser.parse_n_electrons()
        if electrons is None:
            raise RuntimeError(f"number of electrons not found in {output}")
        return _Reference(output, orbitals, parser.parse_energy(), *electrons,
                          states=parser.parse_excited_states(mult), reused_step=step)

    # -- target determinant --------------------------------------------------

    EXCITATION_FILE = "excitation.json"

    @classmethod
    def _load_inherited(cls, conf, step) -> Dict:
        path = Path(conf.folder) / f"protocol_{step}" / "deltascf" / cls.EXCITATION_FILE
        if not path.exists():
            raise RuntimeError(f"excitation 'inherit': step {step} has no accepted ΔSCF state for this conformer")
        data = json.loads(path.read_text())
        if data.get("n") is None:
            raise RuntimeError(f"excitation 'inherit': step {step} used a raw alphaconf/betaconf override")
        return data

    @classmethod
    def _save_excitation(cls, base: Path, exc: Optional[Excitation], gap: float, protocol, energy: float) -> None:
        base.mkdir(parents=True, exist_ok=True)
        data = {
            "spin": exc.spin if exc else None, "n": exc.n if exc else None, "m": exc.m if exc else None,
            "label": exc.label if exc else "raw alphaconf/betaconf", "weight": exc.weight if exc else None,
            "source": exc.source if exc else "raw", "gap_ev": gap, "opt": bool(protocol.opt),
            "energy": energy, "step": int(protocol.number),
        }
        (base / cls.EXCITATION_FILE).write_text(json.dumps(data, indent=2))

    def _compare_inherited_gap(self, gap: float, inherited: Dict, settings, tag: str) -> None:
        src_gap, step = inherited.get("gap_ev"), inherited.get("step")
        self.log.info(f"{tag}: gap {gap:.2f} eV (step {step}: {src_gap:.2f} eV)")
        if inherited.get("opt"):
            return  # adiabatic gap of an optimisation: not comparable
        if abs(gap - src_gap) > settings["inherit_gap_tol_ev"]:
            self.log.warning(
                f"{tag}: gap differs by {gap - src_gap:+.2f} eV from step {step} "
                f"(> {settings['inherit_gap_tol_ev']} eV): the inherited {inherited.get('label')} may "
                "correspond to an orbital of different character at this level"
            )

    @staticmethod
    def _target(settings, ref: _Reference, mult: int, engine,
                inherited: Optional[Dict] = None) -> Tuple[_Target, Optional[Excitation], str]:
        if inherited is not None:
            exc = Excitation(inherited["spin"], int(inherited["n"]), int(inherited["m"]), source="inherit")
            return _Target([exc]), exc, ds.classify(mult, exc)
        if settings.get("alphaconf") or settings.get("betaconf"):
            target = _Target(alphaconf=settings.get("alphaconf"), betaconf=settings.get("betaconf"))
            return target, None, "singlet_ap" if mult == 1 else "raw"
        if settings["excitation"] == "auto":
            exc = ds.resolve_auto(ref.states, settings["root"], ref.n_alpha, ref.n_beta, engine.index_base)
        else:
            n, m = ds.parse_excitation(settings["excitation"])
            exc = Excitation(settings["spin"], n, m)
        return _Target([exc]), exc, ds.classify(mult, exc)

    # -- convergence ladder -------------------------------------------------

    def _ladder(self, engine, settings, protocol, ref, target, base, geometry, kind, tag) -> Optional[_Run]:
        cation = None
        for i, attempt in enumerate(ds.ladder(engine.name, settings), start=1):
            guess = ref.orbitals
            if attempt.guess == "cation":
                if cation is None:
                    cation = self._cation(engine, protocol, ref, base, geometry)
                if not cation.ok or cation.orbitals is None:
                    self.log.warning(f"{tag}: cation guess failed ({cation.reason}), skipping attempt {i}")
                    continue
                guess = cation.orbitals
            calc = engine.build(base / f"try{i}", mode="protocol", guess=guess, target=target,
                                attempt=attempt, n_occ=(ref.n_alpha, ref.n_beta))
            result = engine.run(calc, geometry, opt=bool(protocol.opt))
            reason = result.reason if not result.ok else self._reject(result, ref, settings, protocol.mult, kind)
            if reason:
                self.log.info(f"{tag}: attempt {i} ({attempt.name}) rejected: {reason}")
                continue
            gap = (result.energy - ref.energy) * ds.EH_TO_EV
            self.log.info(
                f"{tag}: attempt {i} ({attempt.name}) accepted: E = {result.energy:.6f} Eh, "
                f"gap {gap:.2f} eV" + (f", <S2> = {result.s2:.3f}" if result.s2 is not None else "")
            )
            self._identity_check(gap, ref, settings, tag)
            return result
        return None

    def _cation(self, engine, protocol, ref, base, geometry) -> _Run:
        calc = engine.build(base / "cation", mode="sp", charge=protocol.charge + 1,
                            mult=ds.cation_multiplicity(protocol.mult), guess=ref.orbitals)
        return engine.run(calc, geometry, opt=False)

    def _reject(self, result: _Run, ref: _Reference, settings, mult: int, kind: str) -> Optional[str]:
        gap = (result.energy - ref.energy) * ds.EH_TO_EV
        if gap < settings["min_gap_ev"]:
            return f"collapse to the reference (gap {gap:.3f} eV < {settings['min_gap_ev']} eV)"
        if kind == "singlet_ap" and result.s2 is not None and result.s2 < 0.2:
            return f"collapse to a closed-shell solution (<S2> = {result.s2:.3f})"
        if result.s2 is not None:
            if kind == "singlet_ap" and not 0.8 <= result.s2 <= 1.2:
                self.log.warning(f"<S2> = {result.s2:.3f} outside [0.8, 1.2] for an open-shell singlet")
            elif kind in ("pure", "raw") and abs(result.s2 - ds.expected_s2(mult)) > 0.1:
                self.log.warning(f"<S2> = {result.s2:.3f}, expected {ds.expected_s2(mult):.2f}: spin contamination")
        return None

    def _identity_check(self, gap: float, ref: _Reference, settings, tag: str) -> None:
        root = settings.get("root", 1)
        energies = {s.root: s.energy_ev for s in ref.states}
        if settings["excitation"] != "auto" or root not in energies:
            return
        lower = [r for r in energies if r < root and abs(energies[root] - energies[r]) > 0.3]
        if lower:
            closest = min(energies, key=lambda r: abs(energies[r] - gap))
            if closest in lower:
                self.log.warning(
                    f"{tag}: gap {gap:.2f} eV is closer to TD-DFT root {closest} ({energies[closest]:.2f} eV) "
                    f"than to root {root} ({energies[root]:.2f} eV): possible collapse onto a lower state"
                )

    # -- spin purification ---------------------------------------------------

    def _purify(self, engine, settings, protocol, ref, exc, kind, accepted, base, geometry, tag) -> float:
        if kind == "singlet_ap":
            return self._purify_singlet(engine, protocol, accepted, base, geometry, tag)
        if kind == "doublet_recoupling":
            return self._recouple(engine, settings, protocol, ref, exc, accepted, base, geometry, tag)
        if kind == "unsupported":
            self.log.warning(
                f"{tag}: {exc.label} opens four or more shells: no validated spin purification "
                f"(<S2> = {accepted.s2}); see Hait & Head-Gordon, arXiv:2103.04573"
            )
        else:
            self.log.info(f"{tag}: spin purification not needed (<S2> = {accepted.s2})")
        return accepted.energy

    def _purify_singlet(self, engine, protocol, accepted, base, geometry, tag) -> float:
        calc = engine.build(base / "triplet", mode="sp", mult=3, guess=accepted.orbitals)
        triplet = engine.run(calc, geometry, opt=False)
        if not triplet.ok or triplet.s2 is None or accepted.s2 is None:
            self.log.warning(f"{tag}: triplet for spin purification failed ({triplet.reason}); mixed energy kept")
            return accepted.energy
        if not 1.9 <= triplet.s2 <= 2.1 or triplet.energy >= accepted.energy:
            self.log.warning(
                f"{tag}: triplet not usable for purification (E_T = {triplet.energy:.6f}, <S2> = {triplet.s2:.3f}); "
                "mixed energy kept"
            )
            return accepted.energy
        e_s = ds.yamaguchi(accepted.energy, accepted.s2, triplet.energy, triplet.s2)
        note = " (at the final geometry of the mixed-state optimisation)" if protocol.opt else ""
        self.log.info(
            f"{tag}: spin-purified singlet E_S = {e_s:.6f} Eh (E_mix = {accepted.energy:.6f}, "
            f"E_T = {triplet.energy:.6f}, {(e_s - accepted.energy) * ds.EH_TO_EV:+.3f} eV){note}"
        )
        return e_s

    def _recouple(self, engine, settings, protocol, ref, exc, accepted, base, geometry, tag) -> float:
        dets, target_name = ds.doublet_determinants(exc)
        energies: Dict[str, float] = {target_name: accepted.energy}
        n_total = ref.n_alpha + ref.n_beta
        for name, (mult, excitations) in dets.items():
            if name == target_name:
                continue
            n_alpha = (n_total + mult - 1) // 2
            result = None
            # every determinant keeps its occupation from the reference orbitals
            # (MOM/IMOM in NWChem also for the quartet: plain aufbau drifts to
            # another, higher quartet)
            for i, attempt in enumerate(ds.ladder(engine.name, {**settings, "cation_guess": False}), start=1):
                calc = engine.build(base / f"recoupling_{name}_try{i}", mode="sp", mult=mult,
                                    guess=ref.orbitals, target=_Target(excitations),
                                    attempt=attempt, n_occ=(n_alpha, n_total - n_alpha))
                result = engine.run(calc, geometry, opt=False)
                if result.ok:
                    break
            if result is None or not result.ok:
                self.log.warning(f"{tag}: recoupling determinant {name} failed; mixed energy kept")
                return accepted.energy
            energies[name] = result.energy
        e_q, e_m1, e_m2, e_m3 = (energies[k] for k in ("Q", "M1", "M2", "M3"))
        couplings = ds.exchange_couplings(e_q, e_m1, e_m2, e_m3)
        if min(couplings) <= 0 or e_q >= min(e_m1, e_m2, e_m3):
            self.log.warning(
                f"{tag}: inconsistent recoupling energies (K = {[round(k, 5) for k in couplings]}); mixed energy kept"
            )
            return accepted.energy
        lower, upper = ds.recouple_doublet(e_q, e_m1, e_m2, e_m3)
        chosen = lower if settings.get("doublet_branch", "lower") == "lower" else upper
        self.log.info(
            f"{tag}: recoupled doublets E = {lower:.6f} / {upper:.6f} Eh "
            f"({(lower - ref.energy) * ds.EH_TO_EV:.2f} / {(upper - ref.energy) * ds.EH_TO_EV:.2f} eV vs D0), "
            f"quartet {e_q:.6f}, mixed {e_m1:.6f} / {e_m2:.6f} / {e_m3:.6f}; using the "
            f"{settings.get('doublet_branch', 'lower')} doublet"
        )
        return chosen

    # -- hand-off to the standard parsing -------------------------------------

    def _finalize(self, conf, protocol, engine, accepted: _Run, energy: float, elapsed: float) -> bool:
        from ensemble_analyzer._parser_parameter import get_conf_parameters

        folder = Path.cwd() / conf.folder / f"protocol_{protocol.number}"
        folder.mkdir(parents=True, exist_ok=True)
        output = folder / f"{conf.number}_p{protocol.number}_{engine.name}.{engine.ext}"
        shutil.copyfile(accepted.output, output)
        if accepted.orbitals is not None:
            shutil.copyfile(accepted.orbitals, folder / (engine.orbital_file or accepted.orbitals.name))

        cfg = self.executor.config
        success = get_conf_parameters(
            conf=conf, number=protocol.number, output=output, p=protocol, time=elapsed,
            temp=cfg.temperature, log=self.log, linear=cfg.linear, cut_off=cfg.cut_off,
            alpha=cfg.alpha, P=cfg.P, energy_override=energy,
        )
        if success:
            data = conf.energies[str(protocol.number)]
            self.log.calculation_success(
                conformer_id=conf.number, protocol_number=protocol.number,
                energy=data.E, gibbs=data.G, frequencies=data.Freq, elapsed_time=elapsed,
            )
        return success
