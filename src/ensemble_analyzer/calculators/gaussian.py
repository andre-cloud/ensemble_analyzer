from ase.calculators.gaussian import Gaussian
from ensemble_analyzer.calculators.base import BaseCalc, register_calculator
from ensemble_analyzer.protocol.dispersion import dispersion_keyword

from typing import Tuple, Any

@register_calculator("gaussian")
class GaussianCalc(BaseCalc):

    def _build_calculator(self) -> Tuple[Any, str]:
        from enan_calculators import get_ase_calculator

        solvent = self.protocol.solvent
        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_gaussian"

        calculator = get_ase_calculator(
            "gaussian",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=self.protocol.functional,
            basis=self.protocol.basis,
            solvent=solvent.solvent if solvent else None,
            smd=bool(solvent and solvent.smd),
            cpu=self.cpu,
            add_input=self.protocol.add_input.strip(),
            label=ase_label,
        )

        disp = dispersion_keyword(self.protocol.dispersion, "gaussian")
        if disp:
            calculator.parameters["extra"] += f" EmpiricalDispersion={disp}"

        chk_path = self._build_path(
            self.conf.folder, f"protocol_{self.protocol.number}", "gaussian.chk"
        )
        calculator.parameters["chk"] = chk_path

        if self.protocol.read_orbitals:
            oldchk_path = self._build_path(
                self.conf.folder, f"protocol_{self.protocol.read_orbitals}", "gaussian.chk"
            )
            calculator.parameters["oldchk"] = oldchk_path
            if "guess" not in self.protocol.add_input.lower():
                calculator.parameters["extra"] += " guess=read"

        return calculator, "gaussian"

    def _add_opt_keywords(self, calc: Gaussian) -> None:
        if self.protocol.ts:
            ts_opt = " opt=(ts,calcfc,noeigentest)"
            ts_opt_cons = " opt=(ts,modredundant,calcfc,noeigentest)"
            add_input = self.protocol.add_input or ""
            if 'oldchk' not in add_input and 'readfc' not in add_input and 'calcfc' not in add_input:
                src = self._find_hessian_source_protocol()
                if src is not None:
                    if self.protocol.read_orbitals:
                        if str(self.protocol.read_orbitals) == str(src):
                            ts_opt = " opt=(readfc,ts,noeigentest)"
                            ts_opt_cons = " opt=(readfc,ts,modredundant,noeigentest)"
                    else:
                        chk_path = self._build_path(
                            self.conf.folder, f"protocol_{src}", "gaussian.chk"
                        )
                        calc.parameters["oldchk"] = chk_path
                        ts_opt = " opt=(readfc,ts,noeigentest)"
                        ts_opt_cons = " opt=(readfc,ts,modredundant,noeigentest)"
            opt_str = ts_opt
        else:
            opt_str = " opt"

        if not self.protocol.constraints:
            calc.parameters["extra"] += opt_str
        else:
            if self.protocol.ts:
                calc.parameters["extra"] += ts_opt_cons
            else:
                calc.parameters["extra"] += " opt=(modredundant)"
            tag_map = {1: "X", 2: "B", 3: "A", 4: "D"}
            lines = []
            for c in self.protocol.constraints:
                tag = tag_map.get(len(c), "X")
                lines.append(f"{tag} {' '.join(map(str, c))} F")
            redundant = "\n".join(lines)

            if calc.parameters.get("addsec"):
                calc.parameters["addsec"] += redundant
            else:
                calc.parameters["addsec"] = redundant

        if self.protocol.freq:
            calc.parameters["extra"] += " freq=(HPModes,vcd)"

    def _add_freq_keywords(self, calc: Gaussian) -> None:
        calc.parameters["extra"] += " freq=(HPModes,vcd)"

    def _add_nmr_keywords(self, calc: Gaussian) -> None:
        settings = self._nmr_settings()
        if not settings:
            return
        from ensemble_analyzer._nmr.inputs import gaussian_keywords

        route, atoms_line = gaussian_keywords(settings, self.conf.atoms)
        calc.parameters["extra"] += f" {route}"
        if atoms_line:
            previous = calc.parameters.get("addsec")
            calc.parameters["addsec"] = f"{previous}\n\n{atoms_line}" if previous else atoms_line

    def _add_tddft_keywords(self, calc: Gaussian) -> None:
        nroots = self.protocol.nroots
        if isinstance(nroots, int) and nroots > 0:
            td_str = f" td=(nstates={nroots}"
            tda = self.protocol.tda
            if isinstance(tda, bool) and tda:
                td_str += ",tda"
            td_str += ")"
            calc.parameters["extra"] += td_str
