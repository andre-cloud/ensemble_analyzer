import os
import shutil
from pathlib import Path
import re
from ase.calculators.orca import ORCA, OrcaProfile

_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks(text: str) -> tuple[str, str]:
    """
    Split the input text into blocks that go before coordinates (pre) 
    and blocks that go after coordinates (post).
    Uses regex to gracefully separate blocks starting with %, ! or * 
    even if newlines are missing.
    """
    if not text.strip():
        return text, ""
    
    # 1. Remove comments
    lines = text.split('\n')
    clean_lines = [line.split('#')[0] for line in lines]
    clean_text = '\n'.join(clean_lines)
    
    pre = []
    post = []
    
    # 2. Split right before every %, ! or *
    parts = re.split(r'(?=[%!*])', clean_text)
    
    for part in parts:
        stripped = part.strip()
        if not stripped:
            continue
            
        first_token = stripped.split()[0].lower()
        
        if first_token in _POST_COORDS_KEYWORDS:
            post.append(stripped)
        else:
            pre.append(stripped)
            
    return "\n".join(pre), "\n".join(post)


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

        def patched_write_input(*args, **kwargs):
            original(*args, **kwargs)
            inp = Path(calculator.directory) / calculator.template.inputname
            with open(inp, "a") as f:
                f.write("\n" + post + "\n")

        calculator.write_inputfiles = patched_write_input

    return calculator
