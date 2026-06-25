from __future__ import annotations

import sys

from textual.app import App

from ensemble_analyzer.cli.protocol_wizard.app_state import ProtocolState
from ensemble_analyzer.cli.protocol_wizard.screens import MainScreen


class ProtocolApp(App):
    TITLE = "Protocol Wizard"
    CSS = ""  # screens handle their own styles

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
