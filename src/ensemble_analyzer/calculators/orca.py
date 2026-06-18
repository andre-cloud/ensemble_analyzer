from ase.calculators.orca import ORCA, OrcaProfile
from ensemble_analyzer.calculators.base import BaseCalc, register_calculator
import shutil
import os

from typing import Tuple, Any

_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _is_end_line(line: str) -> bool:
    s = line.split('#')[0].strip()
    return s == 'end'

def _split_post_blocks(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    lines = text.split('\n')
    pre = []
    post = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        kw_match = next((kw for kw in _POST_COORDS_KEYWORDS if stripped.startswith(kw)), None)
        if kw_match:
            block_end = len(lines)
            for j in range(i + 1, len(lines)):
                if any(lines[j].strip().startswith(kw) for kw in _POST_COORDS_KEYWORDS):
                    block_end = j
                    break
            last_end = None
            for j in range(block_end - 1, i, -1):
                if _is_end_line(lines[j]):
                    last_end = j
                    break
            if last_end is not None:
                post.extend(lines[i:last_end + 1])
                i = last_end + 1
                continue
        pre.append(line)
        i += 1
    return '\n'.join(pre), '\n'.join(post)

VERSION = None

try:
    ORCA_COMMAND = os.getenv("ORCACOMMAND") or shutil.which("orca")
    orca_profile = OrcaProfile(command=ORCA_COMMAND)
    VERSION = int(os.getenv("ORCAVERSION")[0])
except Exception:
    orca_profile = None


@register_calculator("orca")
class OrcaCalc(BaseCalc):

    VERSION = VERSION if VERSION else 0

    def common_str(self) -> Tuple[str, str, str]:
        if self.protocol.solvent:
            if "xtb" in self.protocol.functional.lower():
                solv = f"ALPB({self.protocol.solvent.solvent})"
            elif self.protocol.solvent.solvent.strip():
                solv = f" {self.protocol.solvent}"
            else:
                solv = f" CPCM"
        else:
            solv = ""

        si = f"{self.protocol.functional} {self.protocol.basis} {solv} nopop"

        raw_input = self.protocol.add_input.replace("CONF", str(self.conf.folder))
        pre, post = _split_post_blocks(raw_input)

        ob = (
            f"%pal nprocs {self.cpu} end "
            + pre
            + (" %maxcore 5000" if "maxcore" not in raw_input else "")
        )

        return si, ob, post

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

        if self.constrains:
            tag_map = {1: "C", 2: "B", 3: "A", 4: "D"}
            parts = []
            for c in self.constrains:
                tag = tag_map.get(len(c), "C")
                parts.append(f"{{{tag} {' '.join(map(str, c))} C}}")
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
            calc.parameters["orcablocks"] += (
                f"\n%tddft nroots {nroots} tda {tda_val} end\n"
            )

    def _maybe_add_vcd(self, calc: ORCA) -> None:
        if self.VERSION > 5:
            calc.parameters["orcablocks"] += "\n%freq vcd true end\n"
