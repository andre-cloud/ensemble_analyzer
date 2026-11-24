from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Any

from textual.app import App

from ensemble_analyzer.cli.protocol_wizard.screens import MainScreen


@dataclass
class ProtocolState:
    protocol: dict[str, dict[str, Any]] = field(default_factory=dict)
    filename: str | None = None
    dirty: bool = False


class ProtocolApp(App):
    TITLE = "Protocol Wizard"
    # shared by all screens; each screen adds its own layout
    CSS = """
    Screen {
        align: center middle;
    }
    .box {
        border: solid #3b82f6;
        background: #0f172a;
    }
    """

    def __init__(self):
        self.state = ProtocolState()
        super().__init__()

    def on_mount(self):
        self.push_screen(MainScreen())


def main() -> int:
    app = ProtocolApp()
    app.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
