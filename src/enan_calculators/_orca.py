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
            
    # Unisci i blocchi forzando un a capo dopo ognuno per evitare accavallamenti
    pre_str = "\n".join(pre)
    post_str = "\n".join(post)
    
    return pre_str + "\n" if pre_str else "", post_str + "\n" if post_str else ""


def orca_major_version() -> int | None:
    """Major ORCA version from the ``ORCAVERSION`` environment variable, if set."""
    match = re.match(r"\s*(\d+)", os.getenv("ORCAVERSION") or "")
    return int(match.group(1)) if match else None


def orca_nprocs_group(cpu: int) -> int:
    """Largest proper divisor of ``cpu``: cpu divided by its smallest prime factor.

    8 -> 4, 12 -> 6, 6 -> 3; 1 for cpu <= 1 or a prime cpu.
    """
    if cpu <= 1:
        return 1
    p = 2
    while p * p <= cpu:
        if cpu % p == 0:
            return cpu // p
        p += 1
    return 1


def orca_pal_block(cpu: int, raw_input: str = "", version: int | None = None) -> str:
    """``%pal`` block with ``nprocs`` and, when useful, ``nprocs_group``.

    ``nprocs_group`` lets ORCA >= 5 run independent sub-tasks (numerical
    frequencies, NEB images, ...) in groups. It is omitted when the group
    would be 1, when the user already sets it in ``add_input``, or when the
    ORCA version is known to be older than 5.
    """
    group = orca_nprocs_group(cpu)
    if (
        group > 1
        and "nprocs_group" not in (raw_input or "").lower()
        and (version is None or version >= 5)
    ):
        return f"%pal nprocs {cpu} nprocs_group {group} end\n"
    return f"%pal nprocs {cpu} end\n"


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
        orca_pal_block(cpu, raw_input, orca_major_version())
        + pre
        + ("\n%maxcore 5000\n" if "maxcore" not in raw_input else "")
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
