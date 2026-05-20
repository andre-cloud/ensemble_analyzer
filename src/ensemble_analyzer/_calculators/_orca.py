from ase.calculators.orca import ORCA, OrcaProfile
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator
import shutil
import os
from pathlib import Path

from typing import Tuple

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
    """
    Calculator wrapper for ORCA.
    Handles input generation for SP, OPT, and FREQ jobs.
    """

    label = "orca"
    VERSION = VERSION if VERSION else 0

    def common_str(self) -> Tuple[str, str, str]:
        """
        Generate ORCA simple input, pre-coordinate blocks, and post-coordinate blocks.

        Post-coordinate blocks (%frag, %eprnmr, %nmr, %rel, %epr) are automatically
        extracted from add_input and must appear after the *xyz section in ORCA syntax.

        Returns:
            Tuple[str, str, str]: (simple_input, pre_blocks, post_blocks)
        """

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

        raw_input = self.protocol.add_input.format(CONF=self.conf.folder)
        pre, post = _split_post_blocks(raw_input)

        ob = (
            f"%pal nprocs {self.cpu} end "
            + pre
            + (" %maxcore 5000" if "maxcore" not in raw_input else "")
        )

        return si, ob, post

    def _std_calc(self) -> Tuple[ORCA, str]:
        """
        Create standard ORCA calculator with common settings.

        Post-coordinate blocks (%frag, %eprnmr, etc.) from add_input are
        auto-split and appended after the *xyz section via write_input patching.

        Returns:
            Tuple[ORCA, str]: Initialized ASE ORCA calculator and label.
        """
        si, ob, post = self.common_str()

        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_orca"

        calculator = ORCA(
            profile=orca_profile,
            label=ase_label,
            orcasimpleinput=si,
            orcablocks=ob,
            charge=self.protocol.charge,
            mult=self.protocol.mult,
        )

        if self.protocol.read_orbitals:
            calculator.parameters["orcasimpleinput"] += " moread"
            calculator.parameters[
                "orcablocks"
            ] += f'\n%moinp "{self.conf.folder}/protocol_{self.protocol.read_orbitals}/{self.conf.number}_p{self.protocol.read_orbitals}_orca.gbw"\n'

        if "freq" in self.protocol.add_input.lower():
            calculator.parameters["orcablocks"] += "\n%freq vcd true end\n"

        if post:
            original = calculator.write_input

            def patched_write_input(atoms, properties=None, system_changes=None):
                original(atoms, properties, system_changes)
                inp = Path(calculator.directory) / calculator.input_filename()
                with open(inp, "a") as f:
                    f.write("\n" + post + "\n")

            calculator.write_input = patched_write_input

        return calculator, "orca"

    def single_point(self) -> Tuple[ORCA, str]:
        """Configure Single Point calculation."""
        return self._std_calc()

    def optimisation(self) -> Tuple[ORCA, str]:
        """
        Configure Geometry Optimization.
        Adds constraints if specified in protocol.
        """

        calc, label = self._std_calc()
        calc.parameters["orcasimpleinput"] += " opt"
        if self.constrains:
            tag_map = {1: "C", 2: "B", 3: "A", 4: "D"}
            parts = []
            for c in self.constrains:
                tag = tag_map.get(len(c), "C")
                parts.append(f"{{{tag} {' '.join(map(str, c))} C}}")
            text = "\n%geom Constraints " + " ".join(parts) + " end end\n"
            calc.parameters["orcasimpleinput"] += text
            
        if self.protocol.freq:
            calc.parameters["orcasimpleinput"] += " freq"
            if self.VERSION > 5:
                calc.parameters["orcablocks"] += "\n%freq vcd true end\n"

        return calc, label

    def frequency(self) -> Tuple[ORCA, str]:
        """Configure Frequency calculation."""
        
        calc, label = self._std_calc()
        calc.parameters["orcasimpleinput"] += " freq"
        if self.VERSION > 5:
            calc.parameters["orcablocks"] += "\n%freq vcd true end\n"
        return calc, label