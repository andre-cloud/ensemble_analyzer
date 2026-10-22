"""
Quick-start protocol templates.

System templates live in SYSTEM_TEMPLATES (name → list of step dicts).
User templates are JSON files in ~/.ensemble_analyzer/template/ (or the
directory set by ENAN_TEMPLATES_DIR), one file per template, named after the
file stem. Each file uses the protocol file format ({"0": {...}, "1": {...}})
or a plain list of steps.
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


def load_user_templates(
    directory: Path | None = None,
) -> tuple[dict[str, list[dict[str, Any]]], list[str]]:
    """Return (templates, errors); invalid files are skipped and reported in errors."""
    directory = directory or user_templates_dir()
    templates: dict[str, list[dict[str, Any]]] = {}
    errors: list[str] = []
    if not directory.is_dir():
        return templates, errors
    for path in sorted(directory.glob("*.json")):
        try:
            with open(path, encoding="utf-8") as f:
                templates[path.stem] = _parse_steps(json.load(f))
        except Exception as e:
            errors.append(f"{path.name}: {e}")
    return templates, errors


def apply(steps: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(i): copy.deepcopy(step) for i, step in enumerate(steps)}
