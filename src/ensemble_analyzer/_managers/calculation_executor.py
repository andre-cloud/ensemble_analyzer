from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.protocol import Protocol

from ensemble_analyzer.constants import regex_parsing
from ensemble_analyzer._parser_parameter import get_conf_parameters
from ensemble_analyzer.calculators.base import ML_CALCULATORS
from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants, copy_thermochemical_corrections
from ensemble_analyzer.mode_analysis import NormalModeAnalyzer

import os
import shutil
from pathlib import Path
import time
import numpy as np


class CalculationExecutor:

    def __init__(self, config: CalculationConfig, logger: Logger) -> None:
        self.config = config
        self.logger = logger

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def execute(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        cpu: int | None = None,
    ) -> bool:
        per_job_cpu = cpu if cpu is not None else self.config.cpu
        need_retry = False
        new_geom: np.ndarray | None = None

        for attempt in range(2):
            if attempt > 0:
                conf.energies.data.pop(protocol.number, None)
                conf.last_geometry = new_geom.copy()

            if not self._run_single(
                idx, conf, protocol, per_job_cpu, attempt,
            ):
                return False

            self._log_imaginary_localization(conf, protocol)

            if attempt == 0 and protocol.opt and protocol.freq:
                need_retry, new_geom = self._check_imaginary_and_displace(
                    conf, protocol,
                )
                if need_retry:
                    continue
            break

        if need_retry:
            self.logger.warning(
                f"Conf {conf.number} still has problematic imaginary "
                f"frequencies after displacement – deactivating"
            )
            conf.active = False

        return True

    # ------------------------------------------------------------------
    # Single calculation run (ML or QM)
    # ------------------------------------------------------------------

    def _run_single(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        cpu: int,
        attempt: int,
    ) -> bool:
        if attempt == 0:
            self.logger.calculation_start(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                count=idx,
            )
        else:
            self.logger.info(
                f"  Retry – displaced geometry for Conf {conf.number}"
            )

        is_ml = protocol.calculator.lower() in ML_CALCULATORS

        calc_kwargs = {}
        if is_ml:
            calc_kwargs = dict(
                temperature=self.config.temperature,
                linear=self.config.linear,
                cut_off=self.config.cut_off,
                alpha=self.config.alpha,
                P=self.config.P,
            )

        calc, label = protocol.get_calculator(cpu=cpu, conf=conf, **calc_kwargs)
        atoms = conf.get_ase_atoms(calc)
        (Path(conf.folder) / f"protocol_{protocol.number}").mkdir(parents=True, exist_ok=True)

        os.environ['OMP_NUM_THREADS'] = str(cpu)
        os.environ['MKL_NUM_THREADS'] = str(cpu)
        os.environ['OPENBLAS_NUM_THREADS'] = str(cpu)

        start_time = time.perf_counter()

        energy = None

        with self.logger.track_operation(
            "Single calculation",
            conformer_id=conf.number,
            protocol_number=protocol.number,
        ):
            try:
                if is_ml:
                    if protocol.number not in conf.energies:
                        energy = atoms.get_potential_energy()
                else:
                    try:
                        calc.write_inputfiles(atoms, ['energy'])
                    except AttributeError:
                        calc.write_input(atoms, properties=['energy'])
                    try:
                        calc.template.execute(calc.directory, calc.profile)
                    except AttributeError:
                        calc.execute()
            except Exception as e:
                self.logger.debug(e)
                return False

        elapsed = time.perf_counter() - start_time

        if is_ml:
            return self._finalize_ml(conf, protocol, atoms, energy, elapsed,
                                     attempt)
        return self._finalize_qm(conf, protocol, calc, label, elapsed)

    # ------------------------------------------------------------------
    # ML finalisation
    # ------------------------------------------------------------------

    def _finalize_ml(
        self, conf, protocol, atoms, energy, elapsed, attempt,
    ) -> bool:
        if protocol.number not in conf.energies:
            conf.energies.add(
                protocol.number,
                EnergyRecord(E=energy, time=elapsed),
            )
            compute_rotational_constants(conf, protocol.number)
            m_vec = np.asarray(dipole) if (dipole := atoms.get_dipole_moment()) is not None else np.array([1, 1, 1])
            conf.energies.set(protocol.number, "m_vec", m_vec)
            conf.energies.set(protocol.number, "m",
                              float(np.linalg.norm(m_vec)))
            copy_thermochemical_corrections(conf, protocol.number)
        else:
            elapsed = conf.energies[protocol.number].time or 0

        data = conf.energies[protocol.number]
        if attempt == 0:
            self.logger.calculation_success(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=data.E, gibbs=data.G,
                frequencies=data.Freq,
                elapsed_time=elapsed,
            )
        return True

    # ------------------------------------------------------------------
    # QM finalisation
    # ------------------------------------------------------------------

    def _finalize_qm(self, conf, protocol, calc, label, elapsed) -> bool:
        calc_name = protocol.calculator.lower()
        ext = regex_parsing[calc_name]["ext"]
        output_file = (
            Path.cwd() / conf.folder / f"protocol_{protocol.number}" /
            f'{conf.number}_p{protocol.number}_{label}.{ext}'
        )

        try:
            src = Path(calc.directory).resolve() / calc.template.outputname
        except AttributeError:
            pass
        else:
            if src.exists() and src != output_file:
                shutil.move(src, output_file)

        success = get_conf_parameters(
            conf=conf,
            number=protocol.number,
            output=output_file,
            p=protocol,
            time=elapsed,
            temp=self.config.temperature,
            log=self.logger,
            linear=self.config.linear,
            cut_off=self.config.cut_off,
            alpha=self.config.alpha,
            P=self.config.P,
        )

        if success:
            data = conf.energies[protocol.number]
            self.logger.calculation_success(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=data.E, gibbs=data.G,
                frequencies=data.Freq,
                elapsed_time=elapsed,
            )
        return success

    # ------------------------------------------------------------------
    # Imaginary frequency validation + displacement
    # ------------------------------------------------------------------

    def _check_imaginary_and_displace(
        self,
        conf: Conformer,
        protocol: Protocol,
    ) -> tuple[bool, np.ndarray | None]:
        data = conf.energies[protocol.number]
        freqs = data.Freq
        modes = data.NormalModes

        if not isinstance(freqs, np.ndarray) or freqs.size == 0:
            return False, None
        if not isinstance(modes, np.ndarray) or modes.shape[0] == 0:
            return False, None

        threshold = abs(protocol.neg_freq_threshold)
        neg_idx = np.where(freqs < 0)[0]
        if len(neg_idx) == 0:
            if protocol.ts:
                self.logger.warning(
                    f"Conf {conf.number}: TS with 0 imaginary frequencies – "
                    f"deactivating"
                )
                conf.active = False
            return False, None

        analyzer = NormalModeAnalyzer(
            normal_modes=modes,
            geom=conf.last_geometry,
            atoms=conf.atoms,
        )
        significant, _ = analyzer.classify_negative_freqs(freqs, threshold)

        if protocol.ts and protocol.loc_freq:
            return self._handle_ts_imaginary(
                conf, protocol, freqs, neg_idx, significant, analyzer,
            )
        return self._handle_opt_imaginary(
            conf, protocol, freqs, neg_idx, significant, analyzer,
        )

    def _log_imaginary_localization(
        self,
        conf: Conformer,
        protocol: Protocol,
    ) -> None:
        if protocol.number not in conf.energies:
            return
        data = conf.energies[protocol.number]
        freqs = data.Freq
        modes = data.NormalModes
        if not isinstance(freqs, np.ndarray) or freqs.size == 0:
            return
        if not isinstance(modes, np.ndarray) or modes.shape[0] == 0:
            return
        neg_idx = np.where(freqs < 0)[0]
        if len(neg_idx) == 0:
            return

        threshold = abs(protocol.neg_freq_threshold)
        analyzer = NormalModeAnalyzer(
            normal_modes=modes,
            geom=conf.last_geometry,
            atoms=conf.atoms,
        )
        significant, noise = analyzer.classify_negative_freqs(freqs, threshold)

        if noise:
            self.logger.debug(
                f"Conf {conf.number}: noise imaginary freq(s) ≤ "
                f"{protocol.neg_freq_threshold}i: "
                f"{', '.join(f'{freqs[i]:.2f}' for i in noise)}"
            )
        for i in significant:
            details = analyzer.imag_mode_summary(i, protocol.loc_freq)
            parts = [f"mode {i} ({freqs[i]:.2f})"]
            if "fragments" in details:
                parts.append(
                    ", ".join(f"{k}={v:.1f}%" for k, v in details["fragments"].items())
                )
            if details["top_atoms"]:
                parts.append(
                    "atoms: "
                    +                     ", ".join(f"{idx}:{a}({p:.1f}%)" for idx, a, p in details["top_atoms"])
                )
            self.logger.debug(f"Conf {conf.number}: significant imag " + " | ".join(parts))

    # ------------------------------------------------------------------
    # Case A: opt + freq (no TS)
    # ------------------------------------------------------------------

    def _handle_opt_imaginary(
        self, conf, protocol, freqs, neg_idx, significant, analyzer,
    ) -> tuple[bool, np.ndarray | None]:
        if len(significant) == 0:
            return False, None
        if not protocol.auto_displace:
            self.logger.warning(
                f"Conf {conf.number}: {len(significant)} negative "
                f"freq(s) > {protocol.neg_freq_threshold}i – auto_displace "
                f"disabled, proceeding anyway"
            )
            return False, None

        mode_idx = significant[0]
        new_geom = analyzer.displace_geometry(mode_idx, protocol.displace_scale)
        self.logger.info(
            f"  Displacing along mode {mode_idx} (freq={freqs[mode_idx]:.2f}) "
            f"and re-optimising"
        )
        return True, new_geom

    # ------------------------------------------------------------------
    # Case B: TS
    # ------------------------------------------------------------------

    def _handle_ts_imaginary(
        self, conf, protocol, freqs, neg_idx, significant, analyzer,
    ) -> tuple[bool, np.ndarray | None]:
        n_sig = len(significant)
        if n_sig == 0:
            return False, None

        loc_freq = protocol.loc_freq
        min_ov = protocol.min_overlap
        threshold = protocol.neg_freq_threshold

        sorted_idx = sorted(neg_idx, key=lambda i: abs(freqs[i]), reverse=True)

        def _on_target(mode_i: int) -> bool:
            if not loc_freq:
                return True
            frag = analyzer.localize_mode_fragment(mode_i, loc_freq)
            details = ", ".join(f"{k}={v:.1f}%" for k, v in frag.items())
            self.logger.info(f"  Mode {mode_i} ({freqs[mode_i]:.2f}): {details}")
            return sum(frag.values()) >= min_ov

        largest = sorted_idx[0]
        rank2 = sorted_idx[1] if len(sorted_idx) > 1 else None
        rank2_sig = rank2 is not None and abs(freqs[rank2]) > threshold

        largest_on_target = _on_target(largest)
        rank2_on_target = rank2 is not None and rank2_sig and _on_target(rank2)

        match (n_sig, largest_on_target, rank2_sig, rank2_on_target):
            case (1, True, _, _):
                # B.2: single significant, on target → OK
                return False, None

            case (1, False, _, _):
                # B.3: single significant, NOT on target → deactivate
                top = analyzer.imag_mode_summary(largest, fragments=loc_freq)
                atom_info = "; ".join(f"{s} {a}" for a, s, p in top["top_atoms"])
                self.logger.warning(
                    f"Conf {conf.number}: TS mode NOT on target. "
                    f"Largest displacement on atoms: {atom_info}"
                )
                conf.active = False
                return False, None

            case (_, True, False, _):
                # B.4: largest on target, rank2 is noise/absent → OK
                return False, None

            case (_, True, True, _):
                # B.5: largest on target, rank2 significant → displace rank2
                new_geom = analyzer.displace_geometry(rank2, protocol.displace_scale)
                self.logger.info(
                    f"  TS: mode {largest} on target, displacing along "
                    f"2nd mode {rank2} ({freqs[rank2]:.2f}) and re-optimising"
                )
                return True, new_geom

            case (_, False, _, True):
                # B.6: largest NOT on target, rank2 on target → displace largest
                new_geom = analyzer.displace_geometry(largest, protocol.displace_scale)
                self.logger.info(
                    f"  TS: 2nd mode {rank2} on target, displacing along "
                    f"mode {largest} ({freqs[largest]:.2f}) and re-optimising"
                )
                return True, new_geom

        # Fallback: no TS mode localised on target
        top = analyzer.imag_mode_summary(largest, fragments=loc_freq)
        atom_info = "; ".join(f"{s} {a}" for a, s, p in top["top_atoms"])
        self.logger.warning(
            f"Conf {conf.number}: no TS mode localised on target. "
            f"Largest mode on atoms: {atom_info}"
        )
        conf.active = False
        return False, None
