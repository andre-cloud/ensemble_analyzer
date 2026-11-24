from typing import Tuple, Any

from ase.calculators.orca import ORCA
from enan_calculators import orca_major_version
from ensemble_analyzer.calculators.base import BaseCalc, register_calculator
from ensemble_analyzer.protocol.dispersion import dispersion_keyword


@register_calculator("orca")
class OrcaCalc(BaseCalc):

    VERSION = orca_major_version() or 0

    def _build_calculator(self) -> Tuple[Any, str]:
        from enan_calculators import get_ase_calculator

        if self.protocol.solvent:
            if "xtb" in self.protocol.functional.lower():
                solv_name = f"ALPB({self.protocol.solvent.solvent})"
            elif self.protocol.solvent.solvent.strip():
                solv_name = str(self.protocol.solvent)
            else:
                solv_name = "CPCM"
        else:
            solv_name = None

        raw_input = self.protocol.add_input.replace("CONF", str(self.conf.folder))
        ase_dir = f"{self.conf.folder}/protocol_{self.protocol.number}"

        calculator = get_ase_calculator(
            "orca",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=self.protocol.functional,
            basis=self.protocol.basis,
            solvent=solv_name,
            cpu=self.cpu,
            add_input=raw_input,
            directory=ase_dir,
        )

        disp = dispersion_keyword(self.protocol.dispersion, "orca")
        if disp:
            calculator.parameters["orcasimpleinput"] += f" {disp}"

        if self.protocol.read_orbitals:
            calculator.parameters["orcasimpleinput"] += "\n! moread\n"
            gbw_path = self._build_path(
                self.conf.folder, f"protocol_{self.protocol.read_orbitals}", "orca.gbw"
            )
            calculator.parameters["orcablocks"] += f'\n%moinp "{gbw_path}"\n'

        return calculator, "orca"

    def _add_opt_keywords(self, calc: ORCA) -> None:
        calc.parameters["orcasimpleinput"] += "\n! OptTS \n" if self.protocol.ts else "\n! opt\n"

        if self.protocol.ts:
            add_input = self.protocol.add_input or ""
            if 'inhessname' not in add_input and 'calc_hess' not in add_input:
                src = self._find_hessian_source_protocol()
                if src is not None:
                    hess_path = self._build_path(
                        self.conf.folder, f"protocol_{src}", "orca.hess"
                    )
                    calc.parameters["orcablocks"] += f'\n%geom inhessname "{hess_path}" end\n'
                else:
                    calc.parameters["orcablocks"] += '\n%geom calc_hess true end\n'

        if self.constraints:
            tag_map = {1: "C", 2: "B", 3: "A", 4: "D"}
            parts = []
            for c in self.constraints:
                tag = tag_map.get(len(c), "C")
                parts.append(f"{{{tag} {' '.join(map(str, c))} C}}")
            if parts:  # an empty Constraints block is invalid ORCA input
                text = "\n%geom Constraints " + " ".join(parts) + " end end\n"
                calc.parameters["orcasimpleinput"] += text

        if self.protocol.freq:
            calc.parameters["orcasimpleinput"] += "\n! freq\n"
            self._maybe_add_vcd(calc)

    def _add_freq_keywords(self, calc: ORCA) -> None:
        calc.parameters["orcasimpleinput"] += "\n! freq\n"
        self._maybe_add_vcd(calc)

    def _add_tddft_keywords(self, calc: ORCA) -> None:
        nroots = self.protocol.nroots
        if isinstance(nroots, int) and nroots > 0:
            tda = self.protocol.tda
            tda_val = "true" if (isinstance(tda, bool) and tda) else "false"
            # `iroot 0`: FINAL SINGLE POINT ENERGY is the ground state (E(SCF) +
            # dispersion) instead of the energy of root 1. Not for optimisations
            # (the followed root is the one optimised) or when the user sets iroot.
            add_input = (self.protocol.add_input or "").lower()
            follow_gs = (
                not self.protocol.opt
                and "opt" not in add_input
                and "iroot" not in add_input
            )
            iroot = " iroot 0" if follow_gs else ""
            calc.parameters["orcablocks"] += (
                f"\n%tddft nroots {nroots} tda {tda_val}{iroot} end\n"
            )

    def _maybe_add_vcd(self, calc: ORCA) -> None:
        if self.VERSION > 5:
            calc.parameters["orcablocks"] += "\n%freq vcd true end\n"
