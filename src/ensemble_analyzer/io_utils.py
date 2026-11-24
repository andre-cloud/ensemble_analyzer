from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from typing import Any

import numpy as np


def _serialise(obj: Any) -> Any:
    """Recursively convert non-JSON-serializable objects to plain Python types."""
    if isinstance(obj, np.ndarray):
        return _serialise(obj.tolist())
    if isinstance(obj, np.floating):
        return None if np.isnan(obj) or np.isinf(obj) else float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k) if isinstance(k, (int, np.integer)) else k: _serialise(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialise(item) for item in obj]
    if is_dataclass(obj) and not isinstance(obj, type):
        return _serialise(asdict(obj))
    return obj


def _is_matrix_like(val: list) -> bool:
    """True if *val* is a list of lists (2D array) — formatted compactly."""
    return bool(val) and isinstance(val[0], list)


def _format_value(val: Any, indent: int = 0, step: int = 4) -> str:
    """Recursively format JSON, keeping matrix rows on single lines."""
    prefix = " " * indent
    inner = " " * (indent + step)

    if isinstance(val, dict):
        if not val:
            return "{}"
        items = []
        for k, v in val.items():
            items.append(f"{inner}{json.dumps(k)}: {_format_value(v, indent + step, step)}")
        return "{\n" + ",\n".join(items) + "\n" + prefix + "}"

    if isinstance(val, list):
        if not val:
            return "[]"
        if _is_matrix_like(val):
            rows = ",\n".join(inner + json.dumps(row) for row in val)
            return "[\n" + rows + "\n" + prefix + "]"
        items = []
        for item in val:
            items.append(inner + _format_value(item, indent + step, step))
        return "[\n" + ",\n".join(items) + "\n" + prefix + "]"

    return json.dumps(val)


def write_json(data: Any, fp, indent: int = 4) -> None:
    """Write JSON with compact matrix formatting.

    Matrices (lists of lists) are kept on single inner lines::

        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0]
        ]

    instead of expanding every element onto a separate line.

    Args:
        data: Data to serialize (numpy objects are converted automatically).
        fp: File-like object to write to.
        indent: Spaces per level (default 4).
    """
    fp.write(_format_value(_serialise(data), 0, indent))