from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ProtocolState:
    protocol: dict[str, dict[str, Any]] = field(default_factory=dict)
    filename: str | None = None
    dirty: bool = False
