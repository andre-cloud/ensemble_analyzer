"""
Quick-start protocol templates.

System templates live in SYSTEM_TEMPLATES (name → list of step dicts).
User templates are JSON files in ~/.ensemble_analyzer/template/ (or the
directory set by ENAN_TEMPLATES_DIR), one file per template. The directory
layout is free: sub-folders (also symlinked ones, e.g. a shared group folder)
are scanned recursively and each template is keyed by its path relative to the
templates directory, without extension ("my_opt", "lab/dft/opt_freq").
Each file uses the protocol file format ({"0": {...}, "1": {...}}) or a plain
list of steps.
"""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from typing import Any

SYSTEM_TEMPLATES: dict[str, list[dict[str, Any]]] = {
    "Empty — start from scratch": [],
    "SP (Single Point)": [
        {"calculator": "orca", "functional": "B3LYP", "basis": "def2-SVP"},
    ],
    "OPT (Geometry Optimization)": [
        {"calculator": "orca", "functional": "B3LYP", "basis": "def2-SVP", "opt": True},
    ],
    "OPT+FREQ (Optimization + Frequencies)": [
        {
            "calculator": "orca",
            "functional": "r2SCAN-3c",
            "basis": "def2-mTZVPP",
            "opt": True,
            "freq": True,
        },
    ],
    "SP → OPT+FREQ → SP (full workflow)": [
        {"calculator": "orca", "functional": "b97-3c", "basis": "def2-mTZVP"},
        {
            "calculator": "orca",
            "functional": "r2SCAN-3c",
            "basis": "def2-mTZVPP",
            "opt": True,
            "freq": True,
            "freq_fact": 0.98,
        },
        {"calculator": "orca", "functional": "wB97X-D4", "basis": "def2-QZVPP"},
    ],
    "TS Search (Transition State)": [
        {
            "calculator": "orca",
            "functional": "r2SCAN-3c",
            "basis": "def2-mTZVPP",
            "ts": True,
            "opt": True,
            "freq": True,
        },
    ],
    "TD-DFT (Excited States)": [
        {
            "calculator": "orca",
            "functional": "CAM-B3LYP",
            "basis": "def2-TZVP",
            "nroots": 10,
            "tda": True,
        },
    ],
}


def user_templates_dir() -> Path:
    env_path = os.getenv("ENAN_TEMPLATES_DIR")
    if env_path:
        return Path(env_path).expanduser()
    return Path.home() / ".ensemble_analyzer" / "template"


def _parse_steps(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, dict):
        try:
            items = sorted(data.items(), key=lambda x: int(x[0]))
        except ValueError:
            raise ValueError("step keys must be integers") from None
        steps = [step for _, step in items]
    elif isinstance(data, list):
        steps = data
    else:
        raise ValueError("root must be a JSON object or list")
    if not all(isinstance(step, dict) for step in steps):
        raise ValueError("each step must be a JSON object")
    return steps


def _walk_json(root: Path):
    """Yield .json files under root (top folder first), following symlinks and skipping hidden entries."""
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:  # symlink loop
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for name in sorted(filenames):
            if name.endswith(".json") and not name.startswith("."):
                yield Path(dirpath) / name


def load_user_templates(
    directory: Path | None = None,
) -> tuple[dict[str, list[dict[str, Any]]], list[str]]:
    """Return (templates, errors); invalid files are skipped and reported in errors.

    Keys are paths relative to the directory without extension; templates in the
    top folder come first, then each sub-folder in alphabetical order.
    """
    directory = directory or user_templates_dir()
    templates: dict[str, list[dict[str, Any]]] = {}
    errors: list[str] = []
    if not directory.is_dir():
        return templates, errors
    for path in _walk_json(directory):
        rel = path.relative_to(directory)
        try:
            with open(path, encoding="utf-8") as f:
                templates[rel.with_suffix("").as_posix()] = _parse_steps(json.load(f))
        except Exception as e:
            errors.append(f"{rel.as_posix()}: {e}")
    return templates, errors


def user_template_path(name: str, directory: Path | None = None) -> Path:
    """Path of the user template `name` ("my_opt" or "group/my_opt"), validated to stay inside the directory."""
    directory = directory or user_templates_dir()
    rel = Path(name.strip())
    if not rel.parts or rel.is_absolute() or any(p == ".." or p.startswith(".") for p in rel.parts):
        raise ValueError(f"'{name}' is not a valid template name")
    if rel.suffix != ".json":
        rel = rel.with_name(rel.name + ".json")
    return directory / rel


def save_user_template(
    name: str,
    protocol: dict[str, dict[str, Any]],
    directory: Path | None = None,
    overwrite: bool = False,
) -> Path:
    """Write protocol as a user template, creating sub-folders as needed."""
    path = user_template_path(name, directory)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(protocol, f, indent=4, ensure_ascii=False)
    return path


def apply(steps: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(i): copy.deepcopy(step) for i, step in enumerate(steps)}
