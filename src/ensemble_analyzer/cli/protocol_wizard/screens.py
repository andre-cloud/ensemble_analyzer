from __future__ import annotations

import json
from typing import Any

from textual import on
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen, Screen
from textual.widgets import (
    Button,
    Collapsible,
    Input,
    Label,
    ListItem,
    ListView,
    RichLog,
    Select,
    Static,
)

from ensemble_analyzer.cli.protocol_wizard.fields import (
    DEFAULTS,
    ML_CALCULATORS,
    CALCULATOR_CHOICES,
)
from ensemble_analyzer.cli.protocol_wizard.templates import list_names, apply

# ── helpers ────────────────────────────────────────────────────────

def _step_label(key: str, step: dict) -> str:
    calc = step.get("calculator", "?").upper()
    func = step.get("functional", "?")
    parts = [f"Step {key}: {calc} | {func}"]
    if step.get("opt"):
        parts.append("[OPT]")
    if step.get("freq"):
        parts.append("[FREQ]")
    if step.get("ts"):
        parts.append("[TS]")
    return " ".join(parts)


def _v(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    return str(v)


def _parse_bool(v: str) -> bool:
    return v.strip().lower() in ("yes", "true", "1", "y")


def _clean(step: dict) -> dict:
    cleaned = {}
    for key, value in step.items():
        if value is None:
            continue
        if isinstance(value, str) and not value:
            continue
        if isinstance(value, (list, dict)) and not value:
            continue
        if key == "calculator":
            cleaned[key] = value
        elif key in DEFAULTS and value == DEFAULTS[key]:
            continue
        else:
            cleaned[key] = value
    return cleaned


def _parse_coord_list(raw: str) -> list[list[int]]:
    """Parse '1,2;3,4,5' into [[1,2],[3,4,5]]."""
    result = []
    for group in raw.split(";"):
        group = group.strip()
        if not group:
            continue
        atoms = [int(i.strip()) for i in group.split(",") if i.strip()]
        if 2 <= len(atoms) <= 4:
            result.append(atoms)
    return result


def _fmt_constraints(data: list[list[int]]) -> str:
    return ";".join(",".join(str(a) for a in g) for g in data)


def _fmt_validators(data: list[list]) -> str:
    return "\n".join(f"{v[0]} ||| {v[1]} ||| {v[2]}" for v in data)


def _fmt_solvent(data: dict) -> str:
    return f"{data.get('solvent', '')}  SMD={'yes' if data.get('smd') else 'no'}"


# ── Step editor modal ──────────────────────────────────────────────

class StepEditorScreen(ModalScreen[dict | None]):
    # CSS = """
    # StepEditorScreen {
    #     align: center middle;
    # }
    # #editor-box {
    #     width: 80%;
    #     height: 90%;
    #     border: solid #3b82f6;
    #     background: #0f172a;
    # }
    # #editor-box > VerticalScroll {
    #     scrollbar-gutter: stable;
    #     padding: 1 2;
    # }
    # Label.title {
    #     text-style: bold;
    #     background: #3b82f6;
    #     color: #ffffff;
    #     padding: 0 1;
    #     width: 100%;
    # }
    # Label.field {
    #     margin-top: 1;
    #     text-style: bold;
    # }
    # Select, Input {
    #     width: 100%;
    # }
    # .row {
    #     height: 3;
    # }
    # #editor-buttons {
    #     dock: bottom;
    #     height: 3;
    #     align: center middle;
    # }
    # Collapsible {
    #     margin-top: 1;
    # }
    # .val-row {
    #     height: 3;
    # }
    # .val-remove {
    #     width: 5;
    #     margin-left: 1;
    # }
    # """

    def __init__(self, step: dict, step_num: int, is_new: bool = True):
        self.step = dict(step)
        self.step_num = step_num
        self.is_new = is_new
        self._validators = list(step.get("validators", []))
        super().__init__()

    def compose(self):
        kind = "NEW" if self.is_new else "EDIT"
        with Vertical(id="editor-box"):
            with VerticalScroll():
                yield Label(f"Step {self.step_num} — {kind}", classes="title")

                yield Label("Calculator", classes="field")
                yield Select(
                    CALCULATOR_CHOICES,
                    value=self.step.get("calculator", "orca"),
                    id="f_calculator",
                )

                yield Label("Functional / Method", classes="field")
                yield Input(
                    value=_v(self.step.get("functional")),
                    placeholder="B3LYP",
                    id="f_functional",
                )

                yield Label("Basis set", classes="field")
                yield Input(
                    value=_v(self.step.get("basis")),
                    placeholder="def2-SVP",
                    id="f_basis",
                )

                yield Label("Optimization", classes="field")
                yield Select(
                    [("No", "no"), ("Yes", "yes")],
                    value="yes" if self.step.get("opt") else "no",
                    id="f_opt",
                )

                yield Label("Frequencies", classes="field")
                yield Select(
                    [("No", "no"), ("Yes", "yes")],
                    value="yes" if self.step.get("freq") else "no",
                    id="f_freq",
                )

                yield Label("Transition State (TS)", classes="field")
                yield Select(
                    [("No", "no"), ("Yes", "yes")],
                    value="yes" if self.step.get("ts") else "no",
                    id="f_ts",
                )

                yield Label("Multiplicity", classes="field")
                yield Input(
                    value=str(self.step.get("mult", 1)),
                    id="f_mult",
                )

                yield Label("Charge", classes="field")
                yield Input(
                    value=str(self.step.get("charge", 0)),
                    id="f_charge",
                )

                # ── Solvent ──
                with Collapsible(title="Solvent", collapsed=True):
                    sol = self.step.get("solvent") or {}
                    yield Label("Solvent name", classes="field")
                    yield Input(value=_v(sol.get("solvent")), id="f_solvent_name")
                    yield Label("Use SMD?", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if (isinstance(sol, dict) and sol.get("smd")) else "no",
                        id="f_solvent_smd",
                    )

                # ── TD-DFT ──
                with Collapsible(title="TD-DFT (Excited States)", collapsed=True):
                    yield Label("nroots", classes="field")
                    yield Input(
                        value=_v(self.step.get("nroots")),
                        placeholder="10",
                        id="f_nroots",
                    )
                    yield Label("Tamm-Dancoff (TDA)", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if self.step.get("tda") else "no",
                        id="f_tda",
                    )

                # ── Constraints & Monitoring ──
                with Collapsible(title="Constraints & Monitoring", collapsed=True):
                    yield Label("Constraints (semicolon-separated groups)", classes="field")
                    yield Input(
                        value=_fmt_constraints(self.step.get("constrains", [])),
                        placeholder="1,2; 3,4,5",
                        id="f_constrains",
                    )
                    yield Label("Monitor internals (semicolon-separated groups)", classes="field")
                    yield Input(
                        value=_fmt_constraints(self.step.get("monitor_internals", [])),
                        placeholder="1,2; 3,4,5",
                        id="f_monitor_internals",
                    )

                # ── Pruning & Clustering ──
                with Collapsible(title="Pruning & Clustering", collapsed=True):
                    yield Label("Number of clusters (blank = none)", classes="field")
                    yield Input(
                        value=_v(self.step.get("cluster") if isinstance(self.step.get("cluster"), int) else ""),
                        placeholder="",
                        id="f_cluster",
                    )
                    yield Label("Disable pruning?", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if self.step.get("no_prune") else "no",
                        id="f_no_prune",
                    )

                # ── TS Settings ──
                with Collapsible(title="TS Settings", collapsed=True):
                    yield Label("Localize on atoms (semicolon-separated)", classes="field")
                    yield Input(
                        value=_fmt_constraints(self.step.get("loc_freq", [])),
                        placeholder="1,2,3; 4,5,6",
                        id="f_loc_freq",
                    )
                    yield Label("Min localization %", classes="field")
                    yield Input(
                        value=str(self.step.get("min_localization", 40.0)),
                        id="f_min_localization",
                    )
                    yield Label("Auto-displace", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if self.step.get("auto_displace") else "no",
                        id="f_auto_displace",
                    )
                    yield Label("Displace scale", classes="field")
                    yield Input(
                        value=str(self.step.get("displace_scale", 0.3)),
                        id="f_displace_scale",
                    )
                    yield Label("Neg freq threshold (cm⁻¹)", classes="field")
                    yield Input(
                        value=str(self.step.get("neg_freq_threshold", 20.0)),
                        id="f_neg_freq_threshold",
                    )

                # ── Advanced ──
                with Collapsible(title="Advanced", collapsed=True):
                    yield Label("thrG (kcal/mol)", classes="field")
                    yield Input(value=_v(self.step.get("thrG")), id="f_thrG")
                    yield Label("thrB (cm⁻¹)", classes="field")
                    yield Input(value=_v(self.step.get("thrB")), id="f_thrB")
                    yield Label("thrGMAX (kcal/mol)", classes="field")
                    yield Input(value=_v(self.step.get("thrGMAX")), id="f_thrGMAX")
                    yield Label("Freq scaling factor", classes="field")
                    yield Input(
                        value=str(self.step.get("freq_fact", 1.0)),
                        id="f_freq_fact",
                    )
                    yield Label("Skip failed optimizations?", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if self.step.get("skip_opt_fail") else "no",
                        id="f_skip_opt_fail",
                    )
                    yield Label("Block on retention rate?", classes="field")
                    yield Select(
                        [("No", "no"), ("Yes", "yes")],
                        value="yes" if self.step.get("block_on_retention_rate") else "no",
                        id="f_block_on_retention_rate",
                    )
                    yield Label("Read population from step", classes="field")
                    yield Input(
                        value=_v(self.step.get("read_population")),
                        placeholder="",
                        id="f_read_population",
                    )

                # ── Validators ──
                with Collapsible(title="Output Validators", collapsed=True):
                    yield Label("Current validators:", classes="field")
                    yield Vertical(id="val_list")
                    yield Label("")
                    yield Label("Add new validator:", classes="field")
                    yield Label("Regex (single capture group):")
                    yield Input(
                        placeholder=r"Total:\s*(-?\d+\.\d+)",
                        id="f_val_pattern",
                    )
                    yield Label("Expected value:")
                    yield Input(placeholder="0.0", id="f_val_expected")
                    yield Label("Threshold:")
                    yield Input(placeholder="0.001", id="f_val_threshold")
                    yield Button("Add validator", id="val_add")

                # ── Comment ──
                yield Label("Comment", classes="field")
                yield Input(value=_v(self.step.get("comment")), id="f_comment")

            with Horizontal(id="editor-buttons"):
                yield Button("Save", variant="primary", id="save")
                yield Button("Cancel", id="cancel")

    def _read(self) -> dict:
        s = {}

        s["calculator"] = self.query_one("#f_calculator", Select).value

        func = self.query_one("#f_functional", Input).value.strip()
        if func:
            s["functional"] = func

        basis = self.query_one("#f_basis", Input).value.strip()
        if basis or s.get("calculator") not in ML_CALCULATORS:
            s["basis"] = basis

        s["opt"] = _parse_bool(self.query_one("#f_opt", Select).value)
        s["freq"] = _parse_bool(self.query_one("#f_freq", Select).value)
        s["ts"] = _parse_bool(self.query_one("#f_ts", Select).value)

        try:
            s["mult"] = int(self.query_one("#f_mult", Input).value.strip() or "1")
        except ValueError:
            s["mult"] = 1

        try:
            s["charge"] = int(self.query_one("#f_charge", Input).value.strip() or "0")
        except ValueError:
            s["charge"] = 0

        # solvent
        sol_name = self.query_one("#f_solvent_name", Input).value.strip()
        if sol_name:
            s["solvent"] = {
                "solvent": sol_name,
                "smd": _parse_bool(self.query_one("#f_solvent_smd", Select).value),
            }

        # TD-DFT
        nroots = self.query_one("#f_nroots", Input).value.strip()
        if nroots:
            try:
                s["nroots"] = int(nroots)
            except ValueError:
                pass
        s["tda"] = _parse_bool(self.query_one("#f_tda", Select).value)

        # constraints
        c = _parse_coord_list(self.query_one("#f_constrains", Input).value)
        if c:
            s["constrains"] = c
        m = _parse_coord_list(self.query_one("#f_monitor_internals", Input).value)
        if m:
            s["monitor_internals"] = m

        # clustering
        cluster_raw = self.query_one("#f_cluster", Input).value.strip()
        if cluster_raw:
            try:
                s["cluster"] = int(cluster_raw)
            except ValueError:
                pass
        s["no_prune"] = _parse_bool(self.query_one("#f_no_prune", Select).value)

        # TS
        lf = _parse_coord_list(self.query_one("#f_loc_freq", Input).value)
        if lf:
            s["loc_freq"] = lf
        try:
            s["min_localization"] = float(
                self.query_one("#f_min_localization", Input).value.strip() or "40.0"
            )
        except ValueError:
            pass
        s["auto_displace"] = _parse_bool(
            self.query_one("#f_auto_displace", Select).value
        )
        try:
            s["displace_scale"] = float(
                self.query_one("#f_displace_scale", Input).value.strip() or "0.3"
            )
        except ValueError:
            pass
        try:
            s["neg_freq_threshold"] = float(
                self.query_one("#f_neg_freq_threshold", Input).value.strip() or "20.0"
            )
        except ValueError:
            pass

        # advanced
        for field in ("thrG", "thrB", "thrGMAX"):
            raw = self.query_one(f"#f_{field}", Input).value.strip()
            if raw:
                try:
                    s[field] = float(raw)
                except ValueError:
                    pass
        try:
            s["freq_fact"] = float(
                self.query_one("#f_freq_fact", Input).value.strip() or "1.0"
            )
        except ValueError:
            pass
        s["skip_opt_fail"] = _parse_bool(
            self.query_one("#f_skip_opt_fail", Select).value
        )
        s["block_on_retention_rate"] = _parse_bool(
            self.query_one("#f_block_on_retention_rate", Select).value
        )
        rp = self.query_one("#f_read_population", Input).value.strip()
        if rp:
            s["read_population"] = rp

        # validators
        if self._validators:
            s["validators"] = list(self._validators)

        # comment
        comment = self.query_one("#f_comment", Input).value.strip()
        if comment:
            s["comment"] = comment

        return _clean(s)

    def on_mount(self):
        self._refresh_validators()

    def _refresh_validators(self):
        container = self.query_one("#val_list", Vertical)
        container.remove_children()
        for i, (pat, exp, thr) in enumerate(self._validators):
            container.mount(
                Horizontal(
                    Label(f"  {pat}  →  {exp} ± {thr}  ", id=f"val_label_{i}"),
                    Button("✕", id=f"val_remove_{i}", classes="val-remove"),
                    classes="val-row",
                )
            )

    def on_button_pressed(self, event: Button.Pressed):
        bid = event.button.id or ""
        if bid == "save":
            self.dismiss(self._read())
        elif bid == "cancel":
            self.dismiss(None)
        elif bid == "val_add":
            pat = self.query_one("#f_val_pattern", Input).value.strip()
            exp_raw = self.query_one("#f_val_expected", Input).value.strip()
            thr_raw = self.query_one("#f_val_threshold", Input).value.strip()
            if pat and exp_raw and thr_raw:
                try:
                    self._validators.append([pat, float(exp_raw), float(thr_raw)])
                    self._refresh_validators()
                    self.query_one("#f_val_pattern", Input).value = ""
                    self.query_one("#f_val_expected", Input).value = ""
                    self.query_one("#f_val_threshold", Input).value = ""
                except ValueError:
                    pass  # invalid number, silently ignore
        elif bid.startswith("val_remove_"):
            try:
                idx = int(bid.split("_")[2])
                if 0 <= idx < len(self._validators):
                    self._validators.pop(idx)
                    self._refresh_validators()
            except (IndexError, ValueError):
                pass


# ── Template selector modal ───────────────────────────────────────

class TemplateScreen(ModalScreen[str | None]):
    CSS = """
    TemplateScreen {
        align: center middle;
    }
    #template-box {
        width: 50%;
        height: auto;
        border: solid #3b82f6;
        background: #0f172a;
        padding: 1;
    }
    ListView {
        height: 12;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def compose(self):
        with Vertical(id="template-box"):
            yield Label("Select a template:", classes="title")
            yield ListView(*[ListItem(Static(name)) for name in list_names()], id="template_list")

    def on_list_view_selected(self, event: ListView.Selected):
        if event.list_view.index is not None:
            names = list_names()
            if event.list_view.index < len(names):
                self.dismiss(names[event.list_view.index])


# ── Summary modal ─────────────────────────────────────────────────

class SummaryScreen(ModalScreen[None]):
    CSS = """
    SummaryScreen {
        align: center middle;
    }
    #summary-box {
        width: 80%;
        height: 90%;
        border: solid #3b82f6;
        background: #0f172a;
    }
    #summary-box > VerticalScroll {
        padding: 1 2;
    }
    """

    BINDINGS = [("escape", "dismiss(None)"), ("q", "dismiss(None)")]

    def __init__(self, protocol: dict):
        self.protocol = protocol
        super().__init__()

    def compose(self):
        with Vertical(id="summary-box"):
            with VerticalScroll():
                if not self.protocol:
                    yield Label("(no steps defined)")
                for key in sorted(self.protocol, key=int):
                    s = self.protocol[key]
                    text = _step_label(key, s)
                    yield Label(f"[bold]{text}[/]")
                    for k, v in s.items():
                        if k == "solvent" and isinstance(v, dict):
                            yield Label(f"  solvent: {_fmt_solvent(v)}")
                        elif k in ("constrains", "monitor_internals", "loc_freq") and isinstance(v, list):
                            yield Label(f"  {k}: {_fmt_constraints(v)}")
                        elif k == "validators" and isinstance(v, list):
                            yield Label(f"  validators: {_fmt_validators(v)}")
                        else:
                            yield Label(f"  {k}: {v}")
                    yield Label("")

            yield Button("Close", id="close", variant="primary")

    def on_button_pressed(self, event: Button.Pressed):
        self.dismiss(None)


# ── Confirm dialog modal ──────────────────────────────────────────

class ConfirmScreen(ModalScreen[bool]):
    CSS = """
    ConfirmScreen {
        align: center middle;
    }
    #confirm-box {
        width: 40%;
        height: auto;
        border: solid #f59e0b;
        background: #0f172a;
        padding: 1;
    }
    #confirm-buttons {
        height: auto;
        align: center middle;
        margin-top: 1;
    }
    """

    def __init__(self, message: str):
        self.message = message
        super().__init__()

    def compose(self):
        with Vertical(id="confirm-box"):
            yield Label(self.message)
            with Horizontal(id="confirm-buttons"):
                yield Button("Yes", variant="primary", id="yes")
                yield Button("No", id="no")

    def on_button_pressed(self, event: Button.Pressed):
        self.dismiss(event.button.id == "yes")


# ── Main screen ───────────────────────────────────────────────────

class MainScreen(Screen):
    CSS = """
    MainScreen {
        align: center middle;
    }
    #main-box {
        width: 80%;
        height: 90%;
        border: solid #3b82f6;
        background: #0f172a;
    }
    #header-box {
        height: 3;
        background: #3b82f6;
        padding: 0 1;
        content-align: center middle;
    }
    #body {
        height: 1fr;
    }
    #step-list-box {
        width: 1fr;
        height: 100%;
        border-right: solid #3b82f6;
    }
    #step-list-box > ListView {
        height: 1fr;
    }
    #actions-box {
        width: 28%;
        height: 100%;
        padding: 1;
    }
    #actions-box Button {
        width: 100%;
        margin-bottom: 1;
    }
    #footer-box {
        height: 1;
        background: #0f172a;
        content-align: center middle;
    }
    """

    BINDINGS = [
        Binding("q", "quit_app", "Quit"),
    ]

    def compose(self):
        state = self.app.state  # type: ignore

        with Vertical(id="main-box"):
            yield Static(id="header", classes="header")

            with Horizontal(id="body"):
                with Vertical(id="step-list-box"):
                    yield ListView(id="step_list")

                with Vertical(id="actions-box"):
                    yield Button("New from template", id="new")
                    yield Button("Load protocol", id="load")
                    yield Button("Add step", id="add")
                    yield Button("Edit step", id="edit")
                    yield Button("Remove step", id="remove")
                    yield Button("View summary", id="summary")
                    yield Button("Save", id="save")
                    yield Button("Exit", id="exit")

            yield Static(id="footer", classes="footer")

    def on_mount(self):
        self._refresh()

    def _refresh(self):
        state = self.app.state  # type: ignore
        n = len(state.protocol)
        s = "s" if n != 1 else ""
        unsaved = " ●" if state.dirty else ""
        fname = state.filename or "(new)"
        self.query_one("#header", Static).update(
            f"[bold]Protocol Wizard[/]  |  {fname}  |  {n} step{s}{unsaved}"
        )

        list_view = self.query_one("#step_list", ListView)
        list_view.clear()
        for key in sorted(state.protocol, key=int):
            list_view.append(
                ListItem(Static(_step_label(key, state.protocol[key])))
            )

        self.query_one("#footer", Static).update(
            "[dim]Esc=menu  q=quit[/]"
        )

    def _on_step_edited(self, result: dict | None, key: str | None = None):
        if result is None:
            return
        state = self.app.state  # type: ignore
        if key is not None:
            state.protocol[key] = result
        else:
            k = str(len(state.protocol))
            state.protocol[k] = result
        state.dirty = True
        self._refresh()

    def _on_template_selected(self, name: str | None):
        if name is None:
            return
        state = self.app.state  # type: ignore
        state.protocol = apply(name)
        state.filename = None
        state.dirty = True
        self._refresh()

    # ── Button handlers ──

    def action_quit_app(self):
        state = self.app.state  # type: ignore
        if state.dirty:
            self.app.push_screen(
                ConfirmScreen("Unsaved changes. Discard and quit?"),
                self._on_quit_confirm,
            )
        else:
            self.app.exit()

    def _on_quit_confirm(self, result: bool):
        if result:
            self.app.exit()

    def on_button_pressed(self, event: Button.Pressed):
        state = self.app.state  # type: ignore
        btn = event.button.id

        if btn == "new":
            self.app.push_screen(TemplateScreen(), self._on_template_selected)

        elif btn == "load":
            self._do_load()

        elif btn == "add":
            self.app.push_screen(
                StepEditorScreen(step={}, step_num=len(state.protocol), is_new=True),
                self._on_step_edited,
            )

        elif btn == "edit":
            if not state.protocol:
                return
            self._start_edit_flow()

        elif btn == "remove":
            if not state.protocol:
                return
            self._start_remove_flow()

        elif btn == "summary":
            self.app.push_screen(SummaryScreen(state.protocol))

        elif btn == "save":
            self._do_save_maybe_ask()

        elif btn == "exit":
            if state.dirty:
                self.app.push_screen(
                    ConfirmScreen("Unsaved changes. Save before exiting?"),
                    self._on_exit_confirm,
                )
            else:
                self.app.exit()

    def _on_exit_confirm(self, result: bool):
        if result:
            self._save_and_exit()
        else:
            self.app.exit()

    def _save_and_exit(self):
        state = self.app.state  # type: ignore
        if state.filename:
            self._do_save(state.filename)
            self.app.exit()
        else:
            self.app.push_screen(
                _InputScreen("Save as:", "protocol.json"),
                self._on_save_exit_filename,
            )

    def _on_save_exit_filename(self, fname: str | None):
        if fname:
            self._do_save(fname)
        self.app.exit()

    # ── Edit flow ──

    def _start_edit_flow(self):
        state = self.app.state  # type: ignore
        choices = [(f"{k}: {_step_label(k, state.protocol[k])}", k) for k in sorted(state.protocol, key=int)]
        self.app.push_screen(
            _PickScreen(choices, "Select step to edit:"),
            self._on_edit_key,
        )

    def _on_edit_key(self, key: str | None):
        if key is None:
            return
        state = self.app.state  # type: ignore
        if key not in state.protocol:
            return
        self.app.push_screen(
            StepEditorScreen(
                step=state.protocol[key],
                step_num=int(key),
                is_new=False,
            ),
            lambda r, k=key: self._on_step_edited(r, k),
        )

    # ── Remove flow ──

    def _start_remove_flow(self):
        state = self.app.state  # type: ignore
        choices = [(f"{k}: {_step_label(k, state.protocol[k])}", k) for k in sorted(state.protocol, key=int)]
        self.app.push_screen(
            _PickScreen(choices, "Select step to remove:"),
            self._on_remove_key,
        )

    def _on_remove_key(self, key: str | None):
        if key is None:
            return
        state = self.app.state  # type: ignore
        if key not in state.protocol:
            return
        self.app.push_screen(
            ConfirmScreen(f"Remove step {key}?"),
            lambda confirmed, k=key: self._on_remove_confirm(confirmed, k),
        )

    def _on_remove_confirm(self, confirmed: bool, key: str):
        if not confirmed:
            return
        state = self.app.state  # type: ignore
        del state.protocol[key]
        state.protocol = {
            str(i): state.protocol[k]
            for i, k in enumerate(sorted(state.protocol, key=int))
        }
        state.dirty = True
        self._refresh()

    # ── Load flow ──

    def _do_load(self):
        state = self.app.state  # type: ignore
        if state.dirty:
            self.app.push_screen(
                ConfirmScreen("Unsaved changes. Load and discard?"),
                self._on_load_confirm,
            )
        else:
            self._push_load_dialog()

    def _on_load_confirm(self, result: bool):
        if result:
            self._push_load_dialog()

    def _push_load_dialog(self):
        self.app.push_screen(
            _InputScreen("Load protocol file:", "protocol.json"),
            self._on_load_filename,
        )

    def _on_load_filename(self, fname: str | None):
        if not fname:
            return
        state = self.app.state  # type: ignore
        try:
            with open(fname, encoding="utf-8") as f:
                data = json.load(f)
            if not isinstance(data, dict):
                raise ValueError("Root must be a JSON object")
            steps = {}
            for i, (_, step) in enumerate(sorted(data.items(), key=lambda x: int(x[0]))):
                steps[str(i)] = step
            state.protocol = steps
            state.filename = fname
            state.dirty = False
            self._refresh()
        except Exception as e:
            self.app.push_screen(ConfirmScreen(f"Load error: {e}"))

    # ── Save ──

    def _do_save_maybe_ask(self):
        state = self.app.state  # type: ignore
        if not state.filename:
            self.app.push_screen(
                _InputScreen("Save as:", "protocol.json"),
                self._on_save_filename,
            )
        else:
            self._do_save(state.filename)

    def _on_save_filename(self, fname: str | None):
        if not fname:
            return
        state = self.app.state  # type: ignore
        state.filename = fname
        self._do_save(fname)

    def _do_save(self, fname: str):
        state = self.app.state  # type: ignore
        try:
            with open(fname, "w", encoding="utf-8") as f:
                json.dump(state.protocol, f, indent=4, ensure_ascii=False)
            state.dirty = False
            self._refresh()
        except Exception as e:
            self.app.push_screen(ConfirmScreen(f"Save error: {e}"))


# ── Minimal input screen ──────────────────────────────────────────

class _PickScreen(ModalScreen[str | None]):
    CSS = """
    _PickScreen {
        align: center middle;
    }
    #pick-box {
        width: 50%;
        height: auto;
        border: solid #3b82f6;
        background: #0f172a;
        padding: 1;
    }
    ListView {
        height: 12;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def __init__(self, choices: list[tuple[str, str]], prompt: str = "Select:"):
        self.choices = choices
        self.prompt = prompt
        super().__init__()

    def compose(self):
        with Vertical(id="pick-box"):
            yield Label(self.prompt, classes="title")
            items = [ListItem(Static(label)) for label, _ in self.choices]
            yield ListView(*items, id="pick_list")

    def on_list_view_selected(self, event: ListView.Selected):
        if event.list_view.index is not None:
            idx = event.list_view.index
            if 0 <= idx < len(self.choices):
                self.dismiss(self.choices[idx][1])


class _InputScreen(ModalScreen[str | None]):
    CSS = """
    _InputScreen {
        align: center middle;
    }
    #input-box {
        width: 50%;
        height: auto;
        border: solid #3b82f6;
        background: #0f172a;
        padding: 1;
    }
    #input-buttons {
        height: auto;
        align: center middle;
        margin-top: 1;
    }
    """

    def __init__(self, prompt: str, default: str = ""):
        self.prompt = prompt
        self.default = default
        super().__init__()

    def compose(self):
        with Vertical(id="input-box"):
            yield Label(self.prompt)
            yield Input(value=self.default, id="input_field")
            with Horizontal(id="input-buttons"):
                yield Button("OK", variant="primary", id="ok")
                yield Button("Cancel", id="cancel")

    def on_button_pressed(self, event: Button.Pressed):
        if event.button.id == "ok":
            val = self.query_one("#input_field", Input).value.strip()
            self.dismiss(val if val else None)
        else:
            self.dismiss(None)
