from typing import Tuple, Any
import time
import numpy as np

from .base import BaseCalc


class BaseMlCalc(BaseCalc):
    def common_str(self) -> dict:
        return {}

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        raise NotImplementedError

    def single_point(self) -> Tuple[Any, str]:
        calc = self._get_ml_calculator()
        return calc, self.label

    def _run_vibrations(self, atoms):
        ir_intensities = np.array([])
        try:
            from ase.vibrations import Infrared
            vib = Infrared(atoms, name=f"vib_{self.conf.number}_{self.protocol.number}")
            vib.run()
            ir_intensities = vib.get_ir_intensities()
            freqs = vib.get_frequencies()
        except Exception:
            from ase.vibrations import Vibrations
            vib = Vibrations(atoms, name=f"vib_{self.conf.number}_{self.protocol.number}")
            vib.run()
            freqs = vib.get_frequencies()
        try:
            vib.clean()
        except Exception:
            pass
        return freqs, ir_intensities

    def _compute_thermochemistry(self, energy, scaled_freqs):
        from ensemble_analyzer.rrho import free_gibbs_energy
        rec = self.conf.energies[self.protocol.number]
        pos_freq = scaled_freqs[scaled_freqs > 0]
        if len(pos_freq) > 0 and rec.B_vec is not None:
            try:
                g, zpve, h_val, s_val = free_gibbs_energy(
                    SCF=energy, T=self.temperature, freq=pos_freq,
                    mw=self.conf.weight_mass, B=rec.B_vec,
                    m=self.protocol.mult,
                    linear=self.linear, cut_off=self.cut_off,
                    alpha=self.alpha, P=self.P,
                )
                self.conf.energies.set(self.protocol.number, "G", g)
                self.conf.energies.set(self.protocol.number, "G_E", g - energy)
                self.conf.energies.set(self.protocol.number, "zpve", zpve)
                self.conf.energies.set(self.protocol.number, "H", h_val)
                self.conf.energies.set(self.protocol.number, "S", s_val)
            except Exception:
                pass

    def _store_ir_spectrum(self, raw_freqs, ir_intensities):
        if len(ir_intensities) == 0 or len(raw_freqs) != len(ir_intensities):
            return
        from ensemble_analyzer.conformer.spectral_data import SpectralRecord
        self.conf.graphs_data.add(
            protocol_number=self.protocol.number,
            graph_type="IR",
            record=SpectralRecord(X=raw_freqs, Y=ir_intensities),
        )

    def frequency(self) -> Tuple[Any, str]:
        from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants

        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)

        start = time.perf_counter()
        raw_freqs, ir_intensities = self._run_vibrations(atoms)

        freq_fact = self.protocol.freq_fact
        if freq_fact is None:
            freq_fact = 1.0
        scaled_freqs = raw_freqs * freq_fact

        energy = atoms.get_potential_energy()
        elapsed = time.perf_counter() - start

        self.conf.energies.add(
            self.protocol.number,
            EnergyRecord(E=energy, Freq=scaled_freqs, time=elapsed),
        )

        compute_rotational_constants(self.conf, self.protocol.number)
        self._compute_thermochemistry(energy, scaled_freqs)
        self._store_ir_spectrum(raw_freqs, ir_intensities)

        return calc, self.label

    def optimisation(self) -> Tuple[Any, str]:
        from ase.optimize import BFGS
        from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants

        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)

        start = time.perf_counter()

        with BFGS(atoms) as opt:
            opt.run(fmax=0.05)
        self.conf.last_geometry = atoms.get_positions().copy()

        if self.protocol.freq:
            raw_freqs, ir_intensities = self._run_vibrations(atoms)

            freq_fact = self.protocol.freq_fact
            if freq_fact is None:
                freq_fact = 1.0
            scaled_freqs = raw_freqs * freq_fact

            energy = atoms.get_potential_energy()
            elapsed = time.perf_counter() - start

            self.conf.energies.add(
                self.protocol.number,
                EnergyRecord(E=energy, Freq=scaled_freqs, time=elapsed),
            )

            compute_rotational_constants(self.conf, self.protocol.number)
            self._compute_thermochemistry(energy, scaled_freqs)
            self._store_ir_spectrum(raw_freqs, ir_intensities)

        return calc, self.label
