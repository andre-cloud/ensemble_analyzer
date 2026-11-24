"""NMR reference molecules: optimisation and shielding, cached across steps and restarts.

For every reference of an NMR step:

1. the geometry is optimised at the level of ``nmr['reference_opt_step']`` (the last ground-state
   optimisation before the step) or, without one, at the level of the NMR step itself. A geometry
   already optimised at that step and level is reused; otherwise the optimisation starts from the
   last cached geometry (or from the xyz / bundled / SMILES geometry the first time);
2. its shielding is computed at the level of the NMR step (reused when already done);
3. σ_ref of each nucleus is the mean over the reference atoms of that element (or ``atoms``).

Files: ``nmr_references/<name>/opt/protocol_<k>/``, ``nmr_references/<name>/nmr/protocol_<n>/``
and the cache ``nmr_references/<name>/reference.json``.
"""

import copy
import json
import re
from pathlib import Path
from typing import Dict

import numpy as np

from ensemble_analyzer._managers.calculation_executor import CalculationExecutor
from ensemble_analyzer._nmr.nuclei import element_of
from ensemble_analyzer._nmr.references import resolve_reference_geometry
from ensemble_analyzer.conformer.conformer import Conformer

ROOT = Path("nmr_references")
CACHE = "reference.json"


def _level(p) -> str:
    solvent = str(p.solvent) if getattr(p, "solvent", None) else ""
    fields = (p.calculator, p.functional, p.basis, p.dispersion, solvent, (p.add_input or "").strip())
    return " | ".join(str(f or "").strip().lower() for f in fields)


class _CacheOnly:
    """Executor used by enan_regraph: references must already be in the cache."""

    def execute(self, idx, conf, protocol, cpu=None):
        raise RuntimeError(
            f"NMR reference in {conf.folder} was never computed for protocol {protocol.number}"
        )


def cached_reference_sigmas(protocol, protocols, config, logger) -> Dict[str, Dict]:
    """σ_ref of an NMR step from ``nmr_references/*/reference.json``, without running anything."""
    runner = NMRReferenceRunner(config, logger)
    runner.executor = _CacheOnly()
    return runner.ensure(protocol, protocols)


def _safe(name: str) -> str:
    return re.sub(r"[^\w.+-]", "_", name)


class NMRReferenceRunner:

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.executor = CalculationExecutor(config, logger)

    # ------------------------------------------------------------------

    def ensure(self, protocol, protocols) -> Dict[str, Dict]:
        """σ_ref per nucleus of an NMR step: ``{nucleus: {sigma, reference, delta, n_atoms, opt_step}}``."""
        settings = getattr(protocol, "nmr", None)
        if not settings or not settings.get("shielding"):
            return {}
        by_number = {int(p.number): p for p in protocols}
        out = {}
        for ref in settings["references"]:
            out.update(self._one(ref, protocol, by_number))
        return out

    # ------------------------------------------------------------------

    def _one(self, ref: Dict, protocol, by_number) -> Dict[str, Dict]:
        folder = ROOT / _safe(ref["name"])
        folder.mkdir(parents=True, exist_ok=True)
        cache_file = folder / CACHE
        cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
        cache.setdefault("geometries", {})
        cache.setdefault("shielding", {})

        geom0 = resolve_reference_geometry(ref, base_dir=Path.cwd())
        cache.update(name=ref["name"], source=geom0.source, charge=geom0.charge, mult=geom0.mult)

        k = protocol.nmr.get("reference_opt_step")
        opt_source = by_number.get(int(k)) if k is not None else protocol
        opt_key, opt_level = str(opt_source.number), _level(opt_source)

        geo = cache["geometries"].get(opt_key)
        if geo is None or geo.get("level") != opt_level:
            start = self._start_geometry(cache, geom0)
            self.logger.info(
                f"NMR reference {ref['name']}: optimisation at the level of protocol {opt_key} "
                f"({opt_source.functional}/{opt_source.basis})"
            )
            opt = self._variant(opt_source, geom0.charge, geom0.mult, opt=True)
            conf = self._run(ref, folder / "opt", start, opt)
            geo = {"symbols": list(start[0]), "coords": conf.last_geometry.tolist(), "level": opt_level}
            cache["geometries"][opt_key] = geo
            cache["shielding"] = {
                s: v for s, v in cache["shielding"].items() if v.get("opt_step") != opt_key
            }
            self._save(cache_file, cache)

        n_key, n_level = str(protocol.number), _level(protocol)
        shield = cache["shielding"].get(n_key)
        if shield is None or shield.get("level") != n_level or shield.get("opt_step") != opt_key \
                or set(ref["nuclei"]) - set(shield.get("nuclei", [])):
            sp = self._variant(protocol, geom0.charge, geom0.mult, opt=False)
            sp.nmr = {**protocol.nmr, "nuclei": list(ref["nuclei"]), "coupling": False, "atoms": None,
                      "references": [ref]}
            conf = self._run(ref, folder / "nmr", (geo["symbols"], np.array(geo["coords"])), sp)
            record = conf.nmr_data.get(protocol.number)
            if record is None or not record.shielding:
                raise RuntimeError(f"NMR reference {ref['name']}: no shielding found (protocol {n_key})")
            shield = {"level": n_level, "opt_step": opt_key, "nuclei": list(ref["nuclei"]),
                      "per_atom": {str(i): s for i, s in record.shielding.items()}}
            cache["shielding"][n_key] = shield
            self._save(cache_file, cache)

        per_atom = {int(i): s for i, s in shield["per_atom"].items()}
        out = {}
        for nuc in ref["nuclei"]:
            el = element_of(nuc)
            idx = [i for i, s in enumerate(geo["symbols"]) if s.capitalize() == el]
            if ref.get("atoms") is not None:
                idx = [i for i in idx if i in ref["atoms"]]
            values = [per_atom[i] for i in idx if i in per_atom]
            if not values:
                raise RuntimeError(f"NMR reference {ref['name']}: no {nuc} shielding to average")
            out[nuc] = {"sigma": float(np.mean(values)), "reference": ref["name"],
                        "delta": float(ref.get("delta") or 0.0), "n_atoms": len(values),
                        "opt_step": int(opt_key)}
        return out

    # ------------------------------------------------------------------

    @staticmethod
    def _start_geometry(cache, geom0):
        if cache["geometries"]:
            last = cache["geometries"][max(cache["geometries"], key=int)]
            return last["symbols"], np.array(last["coords"])
        return list(geom0.symbols), geom0.coords

    def _variant(self, source, charge, mult, opt: bool):
        p = copy.copy(source)                 # no __post_init__: keeps Solvent etc.
        p.opt, p.freq, p.ts = opt, False, False
        p.charge, p.mult = int(charge), int(mult)
        p.constraints, p.read_orbitals, p.read_population = [], "", None
        p.nroots, p.deltascf, p.nmr = None, None, None
        p.loc_freq, p.auto_displace = None, False
        p.monitor_internals, p.validators = [], []
        if "CONF" in (p.add_input or ""):
            self.logger.warning(
                f"NMR reference: add_input of protocol {p.number} refers to CONF; "
                "it points to the reference folder"
            )
        return p

    def _run(self, ref, folder: Path, geometry, protocol) -> Conformer:
        symbols, coords = geometry
        conf = Conformer(number=0, geom=np.array(coords, dtype=float), atoms=list(symbols), raw=True)
        conf.folder = str(folder)
        folder.mkdir(parents=True, exist_ok=True)
        ok = self.executor.execute(0, conf, protocol, cpu=self.config.cpu)
        if not ok or not conf.active:
            raise RuntimeError(
                f"NMR reference {ref['name']}: calculation failed in {folder} (protocol {protocol.number})"
            )
        return conf

    @staticmethod
    def _save(path: Path, cache: Dict) -> None:
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(cache, indent=2))
        tmp.replace(path)
