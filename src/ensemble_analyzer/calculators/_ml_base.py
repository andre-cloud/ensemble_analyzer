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
        try:
            from ase.vibrations import Infrared
            ir = Infrared(atoms, name=f"vib_{self.conf.number}_{self.protocol.number}")
            ir.run()
            freqs = ir.get_frequencies()
            ir_intensities = ir.intensities.copy()
        except (AttributeError, NotImplementedError):
            from ase.vibrations import Vibrations
            ir = Vibrations(atoms, name=f"vib_{self.conf.number}_{self.protocol.number}")
            ir.run()
            freqs = ir.get_frequencies()
            ir_intensities = np.zeros(len(freqs))
        try:
            ir.clean()
        except Exception:
            pass
        # ASE returns complex when Hessian has negative eigenvalues.
        # Convert: imaginary freq → negative real, real freq → positive real.
        freqs = np.where(
            np.abs(freqs.imag) > 1e-8,
            -freqs.imag,
            freqs.real,
        )
        return freqs, ir_intensities

    def _compute_thermochemistry(self, energy, scaled_freqs):
        from ensemble_analyzer.conformer.energy_data import compute_thermochemistry
        compute_thermochemistry(
            self.conf, self.protocol.number, energy, scaled_freqs,
            self.temperature, self.linear, self.cut_off, self.alpha, self.P,
            self.protocol.mult,
        )

    def frequency(self) -> Tuple[Any, str]:
        from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants
        from ensemble_analyzer.conformer.spectral_data import SpectralRecord

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

        self.conf.graphs_data.add(
            protocol_number=self.protocol.number,
            graph_type='IR',
            record=SpectralRecord(X=scaled_freqs, Y=ir_intensities),
        )

        compute_rotational_constants(self.conf, self.protocol.number)
        self._compute_thermochemistry(energy, scaled_freqs)

        return calc, self.label

    def optimisation(self) -> Tuple[Any, str]:
        from ase.optimize import BFGS
        from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants
        from ensemble_analyzer.conformer.spectral_data import SpectralRecord

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

            self.conf.graphs_data.add(
                protocol_number=self.protocol.number,
                graph_type='IR',
                record=SpectralRecord(X=scaled_freqs, Y=ir_intensities),
            )

            compute_rotational_constants(self.conf, self.protocol.number)
            self._compute_thermochemistry(energy, scaled_freqs)

        return calc, self.label
