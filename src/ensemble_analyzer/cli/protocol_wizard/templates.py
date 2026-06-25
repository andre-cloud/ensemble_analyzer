"""
Quick-start protocol templates.
Add your own by inserting a new entry in TEMPLATES: name → list of step dicts.
"""

from __future__ import annotations

from typing import Any

TEMPLATES: dict[str, list[dict[str, Any]]] = {
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


def list_names() -> list[str]:
    return list(TEMPLATES.keys())


def apply(name: str) -> dict[str, dict[str, Any]]:
    steps = TEMPLATES.get(name, [])
    return {str(i): step for i, step in enumerate(steps)}
