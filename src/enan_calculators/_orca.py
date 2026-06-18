import os
import shutil
from pathlib import Path
from ase.calculators.orca import ORCA, OrcaProfile

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


def create_orca_calc(charge, mult, method, basis, solvent=None, cpu=1, add_input="", directory=None):
    ORCA_COMMAND = os.getenv("ORCACOMMAND") or shutil.which("orca")
    orca_profile = OrcaProfile(command=ORCA_COMMAND)

    solv = ""
    if solvent:
        solv = f" {solvent}"

    si = f"{method} {basis}{solv} nopop"

    raw_input = add_input or ""
    pre, post = _split_post_blocks(raw_input)

    ob = (
        f"%pal nprocs {cpu} end "
        + pre
        + (" %maxcore 5000" if "maxcore" not in raw_input else "")
    )

    calculator = ORCA(
        profile=orca_profile,
        directory=directory or ".",
        orcasimpleinput=si,
        orcablocks=ob,
        charge=charge,
        mult=mult,
    )

    if post:
        original = calculator.write_inputfiles

        def patched_write_input(atoms, properties=None, system_changes=None):
            original(atoms, properties, system_changes)
            inp = Path(calculator.directory) / calculator.template.inputname
            with open(inp, "a") as f:
                f.write("\n" + post + "\n")

        calculator.write_inputfiles = patched_write_input

    return calculator
