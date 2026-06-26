from typing import Any
import time
from pathlib import Path
import numpy as np
import os

from .base import BaseCalc

from ensemble_analyzer.constants import EV_TO_EH
from ensemble_analyzer.conformer.energy_data import EnergyRecord, compute_rotational_constants, compute_thermochemistry
from ensemble_analyzer.conformer.spectral_data import SpectralRecord
from ase.optimize import LBFGS
from ase.vibrations import Infrared, Vibrations
from ase.constraints import FixAtoms, FixInternals
from sella import Sella


class _CappedSella(Sella):
    """Sella with a hard maxstep cap, like LBFGS."""
    def __init__(self, atoms, maxstep=None, **kwargs):
        super().__init__(atoms, **kwargs)
        self._mstep = maxstep
        if self._mstep is not None:
            self.delta = min(self.delta, self._mstep)

    def step(self):
        super().step()
        if self._mstep is not None:
            self.delta = min(self.delta, self._mstep)


class BaseMlCalc(BaseCalc):
    def common_str(self) -> dict:
        return {}

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        raise NotImplementedError

    def _apply_ase_constraints(self, atoms):
        if not self.constrains:
            return
        fix_atoms = []
        bonds = []
        angles = []
        dihedrals = []
        for c in self.constrains:
            if not isinstance(c, (list, tuple)):
                continue
            idx = [i - 1 for i in c]
            if len(c) == 1:
                fix_atoms.append(idx[0])
            elif len(c) == 2:
                bonds.append([None, idx])
            elif len(c) == 3:
                angles.append([None, idx])
            elif len(c) == 4:
                dihedrals.append([None, idx])
        ase_constraints = []
        if fix_atoms:
            ase_constraints.append(FixAtoms(indices=fix_atoms))
        if bonds or angles or dihedrals:
            ase_constraints.append(FixInternals(
                bonds=bonds if bonds else None,
                angles_deg=angles if angles else None,
                dihedrals_deg=dihedrals if dihedrals else None,
            ))
        atoms.set_constraint(ase_constraints)

    def single_point(self) -> tuple[Any, str]:
        calc = self._get_ml_calculator()
        return calc, self.label

    def _run_vibrations(self, atoms):
        try:
            ir = Infrared(atoms, name=f"{self.conf.folder}/protocol_{self.protocol.number}")
            ir.run()
            freqs = ir.get_frequencies()
            ir_intensities = ir.intensities.copy()
        except (AttributeError, NotImplementedError):
            ir = Vibrations(atoms, name=f"{self.conf.folder}/protocol_{self.protocol.number}")
            ir.run()
            freqs = ir.get_frequencies()
            ir_intensities = np.zeros(len(freqs))

        freqs = np.where(
            np.abs(freqs.imag) > 1e-8,
            -freqs.imag,
            freqs.real,
        )

        masses = atoms.get_masses()
        normal_modes = np.zeros((len(freqs), len(atoms), 3))
        for i in range(len(freqs)):
            mode_cart = ir.get_mode(i)
            normal_modes[i] = mode_cart * np.sqrt(masses[:, None])
        ir.clean()

        return freqs, ir_intensities, normal_modes

    def _compute_thermochemistry(self, energy, scaled_freqs):
        compute_thermochemistry(
            self.conf, self.protocol.number, energy, scaled_freqs,
            self.temperature, self.linear, self.cut_off, self.alpha, self.P,
            self.protocol.mult,
        )

    def _post_optimization_vibrations(self, atoms, start):
        raw_freqs, ir_intensities, normal_modes = self._run_vibrations(atoms)
        freq_fact = self.protocol.freq_fact or 1.0
        scaled_freqs = raw_freqs * freq_fact
        energy = atoms.get_potential_energy() * EV_TO_EH
        elapsed = time.perf_counter() - start
        self.conf.energies.add(
            str(self.protocol.number),
            EnergyRecord(
                E=energy, Freq=scaled_freqs, NormalModes=normal_modes,
                time=elapsed, calculator=self.protocol.calculator,
            ),
        )
        compute_rotational_constants(self.conf, str(self.protocol.number))
        self._compute_thermochemistry(energy, scaled_freqs)
        self.conf.graphs_data.add(
            protocol_number=self.protocol.number,
            graph_type='IR',
            record=SpectralRecord(X=scaled_freqs, Y=ir_intensities),
        )

    def frequency(self) -> tuple[Any, str]:
        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)
        start = time.perf_counter()
        self._post_optimization_vibrations(atoms, start)
        return calc, self.label

    def optimisation(self) -> tuple[Any, str]:
        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)
        self._apply_ase_constraints(atoms)
        start = time.perf_counter()
        logdir = Path(self.conf.folder) / f"protocol_{self.protocol.number}"
        os.mkdir(logdir)

        if self.protocol.ts:
            with _CappedSella(atoms, maxstep=self.protocol.maxstep,
                              logfile=str(logdir / "sella.log")) as opt:
                opt.run(fmax=self.protocol.fmax, steps=self.protocol.maxiter)
        else:
            with LBFGS(atoms, maxstep=self.protocol.maxstep,
                       logfile=str(logdir / "lbfgs.log")) as opt:
                opt.run(fmax=self.protocol.fmax, steps=self.protocol.maxiter)

        self.conf.last_geometry = atoms.get_positions().copy()
        if self.protocol.freq:
            atoms.set_constraint(None)
            self._post_optimization_vibrations(atoms, start)
        return calc, self.label
