"""Text fields of the NMR section of the step editor.

References: ``name[ = smiles | file.xyz]: nuclei; ...`` e.g.
``TMS: 1H 13C; PhF = Fc1ccccc1: 13C; mine = cfcl3.xyz: 19F`` (a bundled name such as TMS needs
neither SMILES nor xyz). Decoupling: ``observed: partners; ...`` e.g. ``13C: 1H``.
"""

import re
from typing import Dict, List, Optional


def _split(text: str) -> List[str]:
    return [p.strip() for p in re.split(r"[;\n]", text or "") if p.strip()]


def _words(text: str) -> List[str]:
    return [w for w in re.split(r"[\s,]+", text.strip()) if w]


def parse_references(text: str) -> List[Dict]:
    refs = []
    for part in _split(text):
        if ":" not in part:
            raise ValueError(f"Reference '{part}': write 'name[=smiles|file.xyz]: nuclei'")
        head, nuclei = part.rsplit(":", 1)
        name, _, source = (x.strip() for x in head.partition("="))
        if not name or not _words(nuclei):
            raise ValueError(f"Reference '{part}': name and nuclei are required")
        ref: Dict = {"name": name}
        if source:
            ref["xyz" if source.lower().endswith(".xyz") else "smiles"] = source
        ref["nuclei"] = _words(nuclei)
        refs.append(ref)
    return refs


def format_references(refs: Optional[List[Dict]]) -> str:
    out = []
    for r in refs or []:
        src = r.get("xyz") or r.get("smiles")
        head = f"{r['name']} = {src}" if src else r["name"]
        out.append(f"{head}: {' '.join(r.get('nuclei', []))}")
    return "; ".join(out)


def parse_decouple(text: str) -> Dict[str, List[str]]:
    out = {}
    for part in _split(text):
        obs, _, partners = part.partition(":")
        if obs.strip() and _words(partners):
            out[obs.strip()] = _words(partners)
    return out


def format_decouple(dec: Optional[Dict[str, List[str]]]) -> str:
    return "; ".join(f"{k}: {' '.join(v)}" for k, v in (dec or {}).items())


def parse_indices(text: str) -> Optional[List[int]]:
    words = _words(text or "")
    return [int(w) for w in words] if words else None


def parse_groups(text: str) -> List[List[int]]:
    """``'1,2,3; 7 8'`` → ``[[1, 2, 3], [7, 8]]``."""
    return [g for g in (parse_indices(p) for p in _split(text)) if g]


def format_groups(groups: Optional[List[List[int]]]) -> str:
    return "; ".join(" ".join(map(str, g)) for g in groups or [])
