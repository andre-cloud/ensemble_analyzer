from ase.calculators.gaussian import Gaussian
from ensemble_analyzer.calculators.base import BaseCalc, register_calculator

from typing import Tuple, Any

@register_calculator("gaussian")
class GaussianCalc(BaseCalc):
    """
    ASE-compatible Gaussian calculator wrapper.
    Encapsulates logic for SP, OPT, and FREQ input generation.
    """

    def common_str(self)-> str:
        """
        Build the common Gaussian route section.

        Returns:
            str: Route string (e.g. '# B3LYP/6-31G* SCRF=...').
        """

        solv = ""
        if self.protocol.solvent:
            if self.protocol.solvent.smd:
                solv = f" SCRF=(SMD,Solvent={self.protocol.solvent.solvent})"
            else:
                solv = f" SCRF=(CPCM,Solvent={self.protocol.solvent.solvent})"

        # Basic route section
        route = f"# {self.protocol.functional}/{self.protocol.basis}{solv}"

        # Add user-specified custom input
        if self.protocol.add_input.strip():
            route += " " + self.protocol.add_input.strip()

        if self.protocol.read_orbitals:
            route += " guess=read"

        return route

    def _build_calculator(self) -> Tuple[Any, str]:
        route = self.common_str()

        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_gaussian"
        label = "gaussian"

        chk_path = self._build_path(
            self.conf.folder, f"protocol_{self.protocol.number}", "gaussian.chk"
        )

        calc = Gaussian(
            label=ase_label,
            output_type='N',
            mem=f"{self.cpu*2}GB",
            chk=chk_path,
            extra=route,
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            nprocshared=self.cpu,
        )
        if self.protocol.read_orbitals:
            oldchk_path = self._build_path(
                self.conf.folder, f"protocol_{self.protocol.read_orbitals}", "gaussian.chk"
            )
            calc.oldchk = oldchk_path

        return calc, label

    def _add_opt_keywords(self, calc: Gaussian) -> None:
        if self.protocol.ts:
            ts_opt = " opt=(ts,calcfc,noeigentest)"
            ts_opt_cons = " opt=(ts,modredudant,calcfc,noeigentest)"
            add_input = self.protocol.add_input or ""
            if 'oldchk' not in add_input and 'readfc' not in add_input and 'calcfc' not in add_input:
                src = self._find_hessian_source_protocol()
                if src is not None:
                    if self.protocol.read_orbitals:
                        if str(self.protocol.read_orbitals) == str(src):
                            ts_opt = " opt=(readfc,ts,noeigentest)"
                            ts_opt_cons = " opt=(readfc,ts,modredudant,noeigentest)"
                    else:
                        chk_path = self._build_path(
                            self.conf.folder, f"protocol_{src}", "gaussian.chk"
                        )
                        calc.oldchk = chk_path
                        ts_opt = " opt=(readfc,ts,noeigentest)"
                        ts_opt_cons = " opt=(readfc,ts,modredudant,noeigentest)"
            opt_str = ts_opt
        else:
            opt_str = " opt"

        if not self.protocol.constrains:
            calc.parameters["extra"] += opt_str
        else:
            if self.protocol.ts:
                calc.parameters["extra"] += ts_opt_cons
            else:
                calc.parameters["extra"] += " opt=(modredudant)"
            tag_map = {1: "X", 2: "B", 3: "A", 4: "D"}
            lines = []
            for c in self.protocol.constrains:
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

    def _add_tddft_keywords(self, calc: Gaussian) -> None:
        nroots = self.protocol.nroots
        if isinstance(nroots, int) and nroots > 0:
            td_str = f" td=(nstates={nroots}"
            tda = self.protocol.tda
            if isinstance(tda, bool) and tda:
                td_str += ",tda"
            td_str += ")"
            calc.parameters["extra"] += td_str
