from __future__ import annotations

import json
from typing import Any

from textual import on
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
    Select,
    Static,
)

from ensemble_analyzer.cli.protocol_wizard.fields import (
    DEFAULTS,
    ML_CALCULATORS,
    CALCULATOR_CHOICES,
)
from ensemble_analyzer.cli.protocol_wizard.nmr_fields import (
    format_decouple, format_groups, format_references, parse_decouple, parse_groups, parse_indices,
    parse_references,
)
from ensemble_analyzer.cli.protocol_wizard.steps import (
    AUTO_KEY,
    apply_defaults,
    as_dict,
    as_list,
    default_read_orbitals,
    export_protocol,
    infer_auto,
    insert_step,
    orbitals_warning,
    remove_step,
    replace_step,
)
from ensemble_analyzer.cli.protocol_wizard.templates import (
    SYSTEM_TEMPLATES,
    load_user_templates,
    save_user_template,
    user_template_path,
    user_templates_dir,
)

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
    if step.get("cpu_per_job"):
        parts.append(f"[{step['cpu_per_job']} CPU/job]")
    if step.get("deltascf"):
        parts.append("[ΔSCF]")
    if step.get("nmr"):
        parts.append("[NMR]")
    if str(step.get("read_orbitals") or "").strip():
        parts.append(f"[orb←{step['read_orbitals']}]")
    return " ".join(parts)


def _v(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    return str(v)


def _parse_bool(v: str) -> bool:
    return v.strip().lower() in ("yes", "true", "1", "y")


def _yn(flag: Any, wid: str, yes_first: bool = False) -> Select:
    """A Yes/No dropdown preset from ``flag``."""
    choices = [("No", "no"), ("Yes", "yes")]
    return Select(choices[::-1] if yes_first else choices,
                  value="yes" if flag else "no", id=wid)


def _clean(step: dict) -> dict:
    cleaned = {}
    for key, value in step.items():
        if value is None:
            continue
        if isinstance(value, str) and not value:
            continue
        if isinstance(value, (list, dict)) and not value:
            continue
        if key in DEFAULTS and value == DEFAULTS[key]:
            continue
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


_CLUSTER_TRUE = ("true", "auto", "yes")
_CLUSTER_FALSE = ("false", "no", "")


def _fmt_cluster(value: Any) -> str:
    """Protocol ``cluster`` value -> wizard text: 'true', 'false' or the number of clusters."""
    if value is True:
        return "true"
    if isinstance(value, int) and not isinstance(value, bool) and value > 1:
        return str(value)
    return "false"


def _parse_cluster(raw: str) -> bool | int:
    """Wizard text -> protocol ``cluster`` value: True (auto k), False (off) or k > 1."""
    raw = raw.strip().lower()
    if raw in _CLUSTER_TRUE:
        return True
    if raw in _CLUSTER_FALSE:
        return False
    try:
        k = int(raw)
    except ValueError:
        return False
    # k < 2 means no clustering (Protocol.clustering requires k > 1)
    return k if k > 1 else False


def _parse_cpu_per_job(raw: str) -> int | None:
    """Wizard text -> protocol ``cpu_per_job`` (None = one job at a time on all CPUs)."""
    try:
        n = int(raw.strip())
    except ValueError:
        return None
    return n if n >= 1 else None


def _fmt_validators(data: list[list]) -> str:
    return "\n".join(f"{v[0]} ||| {v[1]} ||| {v[2]}" for v in data)


def _fmt_solvent(data: dict) -> str:
    return f"{data.get('solvent', '')}  SMD={'yes' if data.get('smd') else 'no'}"


# ── Step editor modal ──────────────────────────────────────────────

class StepEditorScreen(ModalScreen[dict | None]):
    CSS = """
    #editor-box {
        width: 80%;
        height: 90%;
    }
    #editor-box > VerticalScroll {
        scrollbar-gutter: stable;
        padding: 1 2;
    }
    Label.title {
        text-style: bold;
        background: #3b82f6;
        color: #ffffff;
        padding: 0 1;
        width: 100%;
    }
    Label.field {
        margin-top: 1;
        text-style: bold;
    }
    Select, Input {
        width: 100%;
    }
    .row {
        height: 3;
    }
    #editor-buttons {
        dock: bottom;
        height: 3;
        align: center middle;
    }
    Collapsible {
        margin-top: 1;
    }
    .val-row {
        height: 3;
    }
    .val-remove {
        width: 5;
        margin-left: 1;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def __init__(self, step: dict, step_num: int, is_new: bool = True,
                 previous: list[dict] | None = None):
        """``previous``: the steps before this one (for the read_orbitals default)."""
        self.step = dict(step)
        self.step_num = step_num
        self.is_new = is_new
        self.previous = list(previous or [])
        self._validators = list(step.get("validators", []))
        self._ro_default = default_read_orbitals(self.previous, self.step)
        super().__init__()

    def _in(self, wid: str) -> str:
        return self.query_one(f"#{wid}", Input).value.strip()

    def _bool(self, wid: str) -> bool:
        return _parse_bool(self.query_one(f"#{wid}", Select).value)

    def _num(self, wid: str, cast=float, default: str = "", fallback=None):
        """Number typed in an Input: blank -> ``default`` (None if empty);
        not a number -> ``fallback``."""
        raw = self._in(wid) or default
        if not raw:
            return None
        try:
            return cast(raw)
        except ValueError:
            return fallback

    def _put_num(self, s: dict, field: str, default: str = "") -> None:
        """``s[field]`` = float typed in ``#f_<field>``, left unset if blank/invalid."""
        value = self._num(f"f_{field}", float, default)
        if value is not None:
            s[field] = value

    def _ro_default_label(self) -> str:
        default = self._ro_default or "none"
        return (f"Read orbitals from step (default: {default} — previous step of the same QM "
                "program; none for the first step, MLIP/TBLite or a different program)")

    def _current_ro_default(self) -> str:
        """read_orbitals default for the calculator / ΔSCF choice now in the form."""
        probe = {"calculator": self.query_one("#f_calculator", Select).value}
        if self._bool("f_dscf"):
            probe["deltascf"] = {"excitation": "auto"}
        return default_read_orbitals(self.previous, probe)

    @on(Select.Changed, "#f_calculator")
    @on(Select.Changed, "#f_dscf")
    def _update_ro_default(self, event: Select.Changed) -> None:
        new_default = self._current_ro_default()
        field = self.query_one("#f_read_orbitals", Input)
        # a value still equal to the old default follows the new one; a manual value stays
        if field.value.strip() == self._ro_default:
            field.value = new_default
        self._ro_default = new_default
        self.query_one("#l_read_orbitals", Label).update(self._ro_default_label())

    def compose(self):
        kind = "NEW" if self.is_new else "EDIT"
        with Vertical(id="editor-box", classes="box"):
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

                yield Label("Dispersion (orca, gaussian, nwchem)", classes="field")
                yield Input(
                    value=_v(self.step.get("dispersion")),
                    placeholder="empty = none (D2, D3ZERO, D3BJ, D4, PFD)",
                    id="f_dispersion",
                )

                yield Label("ML task / head (uma, fairchem, mace)", classes="field")
                yield Input(
                    value=_v(self.step.get("task")),
                    placeholder="empty = model default (omol for UMA)",
                    id="f_task",
                )

                yield Label("Optimization", classes="field")
                yield _yn(self.step.get("opt"), "f_opt")

                yield Label("Frequencies", classes="field")
                yield _yn(self.step.get("freq"), "f_freq")

                yield Label("Transition State (TS)", classes="field")
                yield _yn(self.step.get("ts"), "f_ts")

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

                yield Label(self._ro_default_label(), classes="field", id="l_read_orbitals")
                ro = self.step.get("read_orbitals")
                yield Input(
                    value=self._ro_default if self.is_new and ro is None else _v(ro),
                    placeholder="empty = no initial guess from a previous step",
                    id="f_read_orbitals",
                )

                # ── Solvent ──
                with Collapsible(title="Solvent", collapsed=True):
                    sol = self.step.get("solvent") or {}
                    yield Label("Solvent name", classes="field")
                    yield Input(value=_v(sol.get("solvent")), id="f_solvent_name")
                    yield Label("Use SMD?", classes="field")
                    yield _yn((isinstance(sol, dict) and sol.get("smd")), "f_solvent_smd")

                # ── TD-DFT ──
                with Collapsible(title="TD-DFT (Excited States)", collapsed=True):
                    yield Label("nroots", classes="field")
                    yield Input(
                        value=_v(self.step.get("nroots")),
                        placeholder="10",
                        id="f_nroots",
                    )
                    yield Label("Tamm-Dancoff (TDA)", classes="field")
                    yield _yn(self.step.get("tda"), "f_tda")

                # ── ΔSCF ──
                dscf = self.step.get("deltascf") or {}
                with Collapsible(title="ΔSCF (Excited State, ORCA/NWChem)", collapsed=not dscf):
                    yield Label("Enable ΔSCF", classes="field")
                    yield _yn(dscf, "f_dscf")
                    yield Label(
                        "Reference TD-DFT step = 'Read orbitals from step' above: it must be a step "
                        "with nroots (GS orbitals, E(GS), state composition). Empty = computed at "
                        "this level", classes="field")
                    yield Label("Excitation: auto (dominant transition of a TD-DFT root), inherit "
                                "(same excitation as the previous ΔSCF step, GS SCF reference, no TD-DFT), "
                                "HOMO->LUMO or HOMO-n->LUMO+m", classes="field")
                    yield Input(value=_v(dscf.get("excitation", "auto")), placeholder="auto", id="f_dscf_exc")
                    yield Label("Root (auto), relative to the step multiplicity: mult 1 -> S1; "
                                "mult 3 -> root 1 = T2 (T1 = plain mult-3 step); mult 2 -> D1", classes="field")
                    yield Input(value=_v(dscf.get("root", 1)), placeholder="1", id="f_dscf_root")
                    yield Label("Spin of the excited electron (explicit excitations)", classes="field")
                    yield Select([("alpha", "alpha"), ("beta", "beta")],
                                 value=dscf.get("spin", "alpha"), id="f_dscf_spin")
                    yield Label("Spin purification: open-shell singlet +1 triplet SP (Yamaguchi); "
                                "doublet with 3 open shells +3 SCF (recoupling); high-spin states: not needed",
                                classes="field")
                    yield Select([("Auto (recommended)", "yes"), ("Off", "no")],
                                 value="yes" if dscf.get("spin_purify", True) else "no", id="f_dscf_purify")
                    yield Label("Doublet kept after recoupling (mult 2)", classes="field")
                    yield Select([("lower", "lower"), ("upper", "upper")],
                                 value=dscf.get("doublet_branch", "lower"), id="f_dscf_branch")
                    yield Label("Cation-orbital guess as last attempts", classes="field")
                    yield _yn(dscf.get("cation_guess", True), "f_dscf_cation", yes_first=True)
                    yield Label("Min. gap vs reference [eV] (lowest state of the same multiplicity)",
                                classes="field")
                    yield Input(value=_v(dscf.get("min_gap_ev", 0.1)), placeholder="0.1", id="f_dscf_gap")
                    yield Label("ORCA alphaconf / betaconf override (advanced, ORCA only)", classes="field")
                    yield Input(value=_v(dscf.get("alphaconf")), placeholder="e.g. 0,1", id="f_dscf_aconf")
                    yield Input(value=_v(dscf.get("betaconf")), placeholder="", id="f_dscf_bconf")

                # ── NMR ──
                nmr = self.step.get("nmr") or {}
                with Collapsible(title="NMR (shielding / coupling, ORCA/Gaussian/NWChem)", collapsed=not nmr):
                    yield Label("Enable NMR (single-point step)", classes="field")
                    yield _yn(nmr, "f_nmr")
                    yield Label("Shielding → shifts (needs references) / spin-spin couplings J", classes="field")
                    yield _yn(nmr.get("shielding", True), "f_nmr_shield", yes_first=True)
                    yield _yn(nmr.get("coupling"), "f_nmr_coupl", yes_first=True)
                    yield Label("Nuclei (isotopes)", classes="field")
                    yield Input(value=" ".join(nmr.get("nuclei") or []), placeholder="1H 13C", id="f_nmr_nuclei")
                    yield Label("Only these atoms (0-based, empty = all)", classes="field")
                    yield Input(value=" ".join(map(str, nmr.get("atoms") or [])), placeholder="", id="f_nmr_atoms")
                    yield Label("References: name[=SMILES|file.xyz]: nuclei; ... "
                                "(TMS, CFCl3, H3PO4, CH3NO2, NH3, BF3OEt2, Me2Se, Me4Sn are bundled)",
                                classes="field")
                    yield Input(value=format_references(nmr.get("references")), placeholder="TMS: 1H 13C",
                                id="f_nmr_refs")
                    yield Label("Extra equivalent atoms (0-based groups, ';'-separated)", classes="field")
                    yield Input(value=format_groups(nmr.get("equivalent")), placeholder="10 11; 14 15",
                                id="f_nmr_equiv")
                    yield Label("Decoupling (observed: partners; ...)", classes="field")
                    yield Input(value=format_decouple(nmr.get("decouple")), placeholder="13C: 1H",
                                id="f_nmr_decouple")
                    yield Label("Spectrometer 1H frequency [MHz] / line width [Hz]", classes="field")
                    yield Input(value=_v(nmr.get("frequency", 400)), placeholder="400", id="f_nmr_freq")
                    yield Input(value=_v(nmr.get("fwhm", 1.0)) if not isinstance(nmr.get("fwhm"), dict) else "1.0",
                                placeholder="1.0", id="f_nmr_fwhm")
                    yield Label("Coupling terms", classes="field")
                    yield Select([("All (FC+SD+PSO+DSO)", "all"), ("Fermi contact only", "fc")],
                                 value=nmr.get("coupling_terms", "all"), id="f_nmr_terms")
                    yield Label("Couplings for a shielding-only step: step number, 'none', empty = auto",
                                classes="field")
                    cf = nmr.get("coupling_from")
                    yield Input(value="none" if cf is False else _v(cf), placeholder="", id="f_nmr_cfrom")
                    yield Label("Level for optimising the references: step number, empty = last opt",
                                classes="field")
                    yield Input(value=_v(nmr.get("opt_from")), placeholder="", id="f_nmr_optfrom")

                # ── Constraints & Monitoring ──
                with Collapsible(title="Constraints & Monitoring", collapsed=True):
                    yield Label("Constraints (semicolon-separated groups)", classes="field")
                    yield Input(
                        value=_fmt_constraints(self.step.get("constraints", [])),
                        placeholder="1,2; 3,4,5",
                        id="f_constraints",
                    )
                    yield Label("Monitor internals (semicolon-separated groups)", classes="field")
                    yield Input(
                        value=_fmt_constraints(self.step.get("monitor_internals", [])),
                        placeholder="1,2; 3,4,5",
                        id="f_monitor_internals",
                    )

                # ── Pruning & Clustering ──
                with Collapsible(title="Pruning & Clustering", collapsed=True):
                    yield Label("Clustering: true (auto k, silhouette), false (off) or number of clusters (needs ≥ 30 conformers)", classes="field")
                    yield Input(
                        value=_fmt_cluster(self.step.get("cluster")),
                        placeholder="true | false | number of clusters",
                        id="f_cluster",
                    )
                    yield Label("Disable pruning?", classes="field")
                    yield _yn(self.step.get("no_prune"), "f_no_prune")
                    yield Label("Duplicate-removal method", classes="field")
                    yield Select(
                        [
                            ("Energy + rotational constants (default)", "rotconst"),
                            ("PRISM: heavy-atom RMSD", "prism"),
                            ("PRISM: RMSD + rotamer-corrected RMSD", "prism_rotcorr"),
                        ],
                        value=self.step.get("prune_method", "rotconst") or "rotconst",
                        id="f_prune_method",
                    )
                    yield Label("PRISM: use moment-of-inertia stage? (opt-in)", classes="field")
                    yield _yn(self.step.get("prism_moi"), "f_prism_moi")
                    yield Label("PRISM: max RMSD (Å)", classes="field")
                    yield Input(
                        value=str(self.step.get("prism_max_rmsd", 0.25)),
                        id="f_prism_max_rmsd",
                    )
                    yield Label("PRISM: max atomic deviation (Å, blank = 2×RMSD)", classes="field")
                    yield Input(
                        value=_v(self.step.get("prism_max_dev")),
                        id="f_prism_max_dev",
                    )
                    yield Label("PRISM: MOI relative deviation (0.01 = 1 %)", classes="field")
                    yield Input(
                        value=str(self.step.get("prism_moi_dev", 0.01)),
                        id="f_prism_moi_dev",
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
                    yield _yn(self.step.get("auto_displace"), "f_auto_displace")
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

                # ── Parallel Execution ──
                with Collapsible(title="Parallel Execution", collapsed=True):
                    yield Label(
                        "CPUs per QM job (default: 1 job at a time on all -cpu cores; "
                        "e.g. 8 with -cpu 32 → 4 jobs at once; ignored for ML)",
                        classes="field",
                    )
                    yield Input(
                        value=_v(self.step.get("cpu_per_job")),
                        placeholder="default: 1 job at a time",
                        id="f_cpu_per_job",
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
                    yield _yn(self.step.get("skip_opt_fail"), "f_skip_opt_fail")
                    yield Label("Block on retention rate?", classes="field")
                    yield _yn(self.step.get("block_on_retention_rate"), "f_block_on_retention_rate")
                    yield Label("Read population from step", classes="field")
                    yield Input(
                        value=_v(self.step.get("read_population")),
                        placeholder="",
                        id="f_read_population",
                    )
                    yield Label("Additional input (raw string)", classes="field")
                    yield Input(
                        value=_v(self.step.get("add_input")),
                        placeholder="int=ultrafine",
                        id="f_add_input",
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

        func = self._in("f_functional")
        if func:
            s["functional"] = func

        basis = self._in("f_basis")
        if basis or s.get("calculator") not in ML_CALCULATORS:
            s["basis"] = basis

        dispersion = self._in("f_dispersion")
        if dispersion:
            s["dispersion"] = dispersion

        task = self._in("f_task")
        if task:
            s["task"] = task

        s["opt"] = self._bool("f_opt")
        s["freq"] = self._bool("f_freq")
        s["ts"] = self._bool("f_ts")

        s["mult"] = self._num("f_mult", int, "1", 1)
        s["charge"] = self._num("f_charge", int, "0", 0)

        # solvent
        sol_name = self._in("f_solvent_name")
        if sol_name:
            s["solvent"] = {
                "solvent": sol_name,
                "smd": self._bool("f_solvent_smd"),
            }

        # TD-DFT
        nroots = self._num("f_nroots", int)
        if nroots is not None:
            s["nroots"] = nroots
        s["tda"] = self._bool("f_tda")

        # ΔSCF
        ref = self._in("f_read_orbitals")
        if self._bool("f_dscf"):
            dscf: dict = {
                "excitation": self._in("f_dscf_exc") or "auto",
                "spin": self.query_one("#f_dscf_spin", Select).value,
                "spin_purify": self._bool("f_dscf_purify"),
                "doublet_branch": self.query_one("#f_dscf_branch", Select).value,
                "cation_guess": self._bool("f_dscf_cation"),
            }
            dscf["root"] = self._num("f_dscf_root", int, "1", 1)
            dscf["min_gap_ev"] = self._num("f_dscf_gap", float, "0.1", 0.1)
            for key, wid in (("alphaconf", "#f_dscf_aconf"), ("betaconf", "#f_dscf_bconf")):
                raw = self.query_one(wid, Input).value.strip()
                if raw:
                    dscf[key] = raw
            s["deltascf"] = dscf
        if ref:
            s["read_orbitals"] = ref

        # NMR
        if self._bool("f_nmr"):
            nmr: dict = {
                "shielding": self._bool("f_nmr_shield"),
                "coupling": self._bool("f_nmr_coupl"),
                "nuclei": self.query_one("#f_nmr_nuclei", Input).value.replace(",", " ").split(),
                "coupling_terms": self.query_one("#f_nmr_terms", Select).value,
            }
            self._nmr_error = None
            try:
                atoms = parse_indices(self.query_one("#f_nmr_atoms", Input).value)
                if atoms is not None:
                    nmr["atoms"] = atoms
                refs = parse_references(self.query_one("#f_nmr_refs", Input).value)
                if refs:
                    nmr["references"] = refs
                groups = parse_groups(self.query_one("#f_nmr_equiv", Input).value)
                if groups:
                    nmr["equivalent"] = groups
                dec = parse_decouple(self.query_one("#f_nmr_decouple", Input).value)
                if dec:
                    nmr["decouple"] = dec
                nmr["frequency"] = float(self._in("f_nmr_freq") or "400")
                nmr["fwhm"] = float(self._in("f_nmr_fwhm") or "1.0")
                cf = self._in("f_nmr_cfrom").lower()
                if cf in ("none", "no", "false"):
                    nmr["coupling_from"] = False
                elif cf:
                    nmr["coupling_from"] = int(cf)
                of = self._in("f_nmr_optfrom")
                if of:
                    nmr["opt_from"] = int(of)
            except ValueError as err:
                self._nmr_error = str(err)
            s["nmr"] = nmr

        # constraints
        c = _parse_coord_list(self.query_one("#f_constraints", Input).value)
        if c:
            s["constraints"] = c
        m = _parse_coord_list(self.query_one("#f_monitor_internals", Input).value)
        if m:
            s["monitor_internals"] = m

        # clustering
        s["cluster"] = _parse_cluster(self.query_one("#f_cluster", Input).value)
        s["no_prune"] = self._bool("f_no_prune")
        s["prune_method"] = self.query_one("#f_prune_method", Select).value or "rotconst"
        if s["prune_method"] != "rotconst":
            s["prism_moi"] = self._bool("f_prism_moi")
            for field in ("prism_max_rmsd", "prism_max_dev", "prism_moi_dev"):
                self._put_num(s, field)

        # TS
        lf = _parse_coord_list(self.query_one("#f_loc_freq", Input).value)
        if lf:
            s["loc_freq"] = lf
        self._put_num(s, "min_localization", "40.0")
        s["auto_displace"] = self._bool("f_auto_displace")
        self._put_num(s, "displace_scale", "0.3")
        self._put_num(s, "neg_freq_threshold", "20.0")

        # parallel execution
        cpu_per_job = _parse_cpu_per_job(self.query_one("#f_cpu_per_job", Input).value)
        if cpu_per_job is not None:
            s["cpu_per_job"] = cpu_per_job

        # advanced
        for field in ("thrG", "thrB", "thrGMAX"):
            self._put_num(s, field)
        self._put_num(s, "freq_fact", "1.0")
        s["skip_opt_fail"] = self._bool("f_skip_opt_fail")
        s["block_on_retention_rate"] = self._bool("f_block_on_retention_rate")
        rp = self._in("f_read_population")
        if rp:
            s["read_population"] = rp
        ai = self._in("f_add_input")
        if ai:
            s["add_input"] = ai

        # validators
        if self._validators:
            s["validators"] = list(self._validators)

        # comment
        comment = self._in("f_comment")
        if comment:
            s["comment"] = comment

        s = _clean(s)
        s[AUTO_KEY] = ref == self._ro_default
        return s

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
            step = self._read()
            ref = str(step.get("read_orbitals") or "")
            if ref and (not ref.isdigit() or int(ref) >= self.step_num):
                self.notify(f"Read orbitals: '{ref}' must be the number of an earlier step "
                            f"(0..{self.step_num - 1})" if self.step_num else
                            "Read orbitals: the first step has no earlier step to read from",
                            severity="error", timeout=8)
                return
            warning = orbitals_warning(self.previous, step, ref)
            if warning:
                self.notify(f"Read orbitals: {warning}", severity="warning", timeout=8)
            if step.get("nmr"):
                from ensemble_analyzer.protocol.nmr import normalize_settings as normalize_nmr
                try:
                    if getattr(self, "_nmr_error", None):
                        raise ValueError(f"NMR: {self._nmr_error}")
                    normalize_nmr(step["nmr"], "-", step.get("calculator", "orca"), step.get("add_input", ""),
                                  step.get("nroots"), opt=bool(step.get("opt") or step.get("ts")),
                                  freq=bool(step.get("freq")), deltascf=step.get("deltascf"))
                except ValueError as err:
                    self.notify(str(err), severity="error", timeout=8)
                    return
            if step.get("deltascf"):
                from ensemble_analyzer.protocol.deltascf import normalize_settings
                try:
                    normalize_settings(step["deltascf"], "-", step.get("calculator", "orca"),
                                       step.get("add_input", ""), step.get("nroots"))
                except ValueError as err:
                    self.notify(str(err), severity="error", timeout=8)
                    return
            self.dismiss(step)
        elif bid == "cancel":
            self.dismiss(None)
        elif bid == "val_add":
            pat = self._in("f_val_pattern")
            exp_raw = self._in("f_val_expected")
            thr_raw = self._in("f_val_threshold")
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

# ── Minimal input screen ──────────────────────────────────────────

class _PickScreen(ModalScreen[str | None]):
    CSS = """
    #pick-box {
        width: 50%;
        height: auto;
        padding: 1;
    }
    #pick_list {
        height: auto;
        max-height: 20;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def __init__(self, choices: list[tuple[str, Any]], prompt: str = "Select:"):
        """``choices``: (label, value) rows; rows with value None are disabled headers."""
        self.choices = choices
        self.prompt = prompt
        super().__init__()

    def compose(self):
        with Vertical(id="pick-box", classes="box"):
            yield Label(self.prompt, classes="title")
            items = [ListItem(Static(label), disabled=value is None)
                     for label, value in self.choices]
            yield ListView(*items, id="pick_list")

    def on_list_view_selected(self, event: ListView.Selected):
        idx = event.list_view.index
        if idx is not None and 0 <= idx < len(self.choices) and self.choices[idx][1] is not None:
            self.dismiss(self.choices[idx][1])


class _InputScreen(ModalScreen[str | None]):
    CSS = """
    #input-box {
        width: 50%;
        height: auto;
        padding: 1;
    }
    #input-buttons {
        height: auto;
        align: center middle;
        margin-top: 1;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def __init__(self, prompt: str, default: str = ""):
        self.prompt = prompt
        self.default = default
        super().__init__()

    def compose(self):
        with Vertical(id="input-box", classes="box"):
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


class TemplateScreen(_PickScreen):
    """Template picker: system templates, then user templates grouped by sub-folder."""

    def __init__(self):
        choices: list[tuple[str, list | None]] = [("[bold]System Templates[/]", None)]
        choices += [(f"  {name}", steps) for name, steps in SYSTEM_TEMPLATES.items()]
        user_templates, errors = load_user_templates()
        choices.append((f"[bold]User Templates[/]  [dim]{user_templates_dir()}[/]", None))
        folder = ""
        for key, steps in user_templates.items():
            parent, _, name = key.rpartition("/")
            if parent != folder:
                folder = parent
                choices.append((f"  [bold]{parent}/[/]", None))
            choices.append((f"    {name}" if parent else f"  {name}", steps))
        if not user_templates:
            choices.append(("  [dim](none — add .json files there or use 'Save as template')[/]", None))
        choices += [(f"  [yellow]⚠ {err}[/]", None) for err in errors]
        super().__init__(choices, "Select a template:")


# ── Summary modal ─────────────────────────────────────────────────

class SummaryScreen(ModalScreen[None]):
    CSS = """
    #summary-box {
        width: 80%;
        height: 90%;
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
        with Vertical(id="summary-box", classes="box"):
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
                        elif k in ("constraints", "monitor_internals", "loc_freq") and isinstance(v, list):
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
    #confirm-box {
        width: 40%;
        height: auto;
        border: solid #f59e0b;
        padding: 1;
    }
    #confirm-buttons {
        height: auto;
        align: center middle;
        margin-top: 1;
    }
    """

    BINDINGS = [("escape", "dismiss(None)")]

    def __init__(self, message: str):
        self.message = message
        super().__init__()

    def compose(self):
        with Vertical(id="confirm-box", classes="box"):
            yield Label(self.message)
            with Horizontal(id="confirm-buttons"):
                yield Button("Yes", variant="primary", id="yes")
                yield Button("No", id="no")

    def on_button_pressed(self, event: Button.Pressed):
        self.dismiss(event.button.id == "yes")


# ── Main screen ───────────────────────────────────────────────────

class MainScreen(Screen):
    CSS = """
    #main-box {
        width: 80%;
        height: 90%;
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
        with Vertical(id="main-box", classes="box"):
            yield Static(id="header", classes="header")

            with Horizontal(id="body"):
                with Vertical(id="step-list-box"):
                    yield ListView(id="step_list")

                with VerticalScroll(id="actions-box"):
                    yield Button("New from template", id="new")
                    yield Button("Load protocol", id="load")
                    yield Button("Add step", id="add")
                    yield Button("Insert step", id="insert")
                    yield Button("Edit step", id="edit")
                    yield Button("Remove step", id="remove")
                    yield Button("View summary", id="summary")
                    yield Button("Save", id="save")
                    yield Button("Save as template", id="save_template")
                    yield Button("Exit", id="exit")

            yield Static("[dim]Tab/navigate  Enter=select  Esc=back  q=quit[/]",
                         id="footer", classes="footer")

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

    def _on_step_edited(self, result: dict | None, key: str | None = None):
        if result is None:
            return
        state = self.app.state  # type: ignore
        state.protocol = replace_step(state.protocol, int(key), result)
        state.dirty = True
        self._refresh()

    def _on_step_inserted(self, result: dict | None, idx: int):
        if result is None:
            return
        state = self.app.state  # type: ignore
        state.protocol = insert_step(state.protocol, idx, result)
        state.dirty = True
        self._refresh()

    def _open_new_step(self, idx: int):
        state = self.app.state  # type: ignore
        self.app.push_screen(
            StepEditorScreen(step={}, step_num=idx, is_new=True,
                             previous=as_list(state.protocol)[:idx]),
            lambda r, i=idx: self._on_step_inserted(r, i),
        )

    def _on_template_selected(self, steps: list | None):
        if steps is None:
            return
        state = self.app.state  # type: ignore
        state.protocol = apply_defaults(as_dict(steps))
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
            self._open_new_step(len(state.protocol))

        elif btn == "insert":
            if not state.protocol:
                self._open_new_step(0)
                return
            choices = [(f"Before {_step_label(k, state.protocol[k])}", k)
                       for k in sorted(state.protocol, key=int)]
            self.app.push_screen(
                _PickScreen(choices, "Insert the new step before:"),
                lambda k: self._open_new_step(int(k)) if k is not None else None,
            )

        elif btn == "edit" and state.protocol:
            self._pick_step("Select step to edit:", self._on_edit_key)

        elif btn == "remove" and state.protocol:
            self._pick_step("Select step to remove:", self._on_remove_key)

        elif btn == "summary":
            self.app.push_screen(SummaryScreen(export_protocol(state.protocol)))

        elif btn == "save":
            self._save()

        elif btn == "save_template":
            if not state.protocol:
                self.app.push_screen(ConfirmScreen("No steps to save as template."))
                return
            self.app.push_screen(
                _InputScreen(
                    f"Template name, sub-folders allowed (e.g. my_opt, group/my_opt)\nSaved in {user_templates_dir()}",
                    "",
                ),
                self._on_template_name,
            )

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
            self._save(then=self.app.exit)
        else:
            self.app.exit()

    def _pick_step(self, prompt: str, callback) -> None:
        state = self.app.state  # type: ignore
        choices = [(f"{k}: {_step_label(k, state.protocol[k])}", k) for k in sorted(state.protocol, key=int)]
        self.app.push_screen(_PickScreen(choices, prompt), callback)

    # ── Edit flow ──

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
                previous=as_list(state.protocol)[:int(key)],
            ),
            lambda r, k=key: self._on_step_edited(r, k),
        )

    # ── Remove flow ──

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
        state.protocol, lost = remove_step(state.protocol, int(key))
        for msg in lost:
            self.notify(msg, severity="warning", timeout=8)
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
            state.protocol = infer_auto(as_dict(as_list(data)))
            state.filename = fname
            state.dirty = False
            self._refresh()
        except Exception as e:
            self.app.push_screen(ConfirmScreen(f"Load error: {e}"))

    # ── Save ──

    def _save(self, then=None):
        """Save to the current file, asking for a name first if there is none;
        ``then`` runs afterwards (also when the name prompt is cancelled)."""
        state = self.app.state  # type: ignore

        def done(fname: str | None):
            if fname:
                state.filename = fname
                self._do_save(fname)
            if then:
                then()

        if state.filename:
            done(state.filename)
        else:
            self.app.push_screen(_InputScreen("Save as:", "protocol.json"), done)

    def _do_save(self, fname: str):
        state = self.app.state  # type: ignore
        try:
            with open(fname, "w", encoding="utf-8") as f:
                json.dump(export_protocol(state.protocol), f, indent=4, ensure_ascii=False)
            state.dirty = False
            self._refresh()
        except Exception as e:
            self.app.push_screen(ConfirmScreen(f"Save error: {e}"))

    # ── Save as user template ──

    def _on_template_name(self, name: str | None):
        if not name:
            return
        try:
            path = user_template_path(name)
        except ValueError as e:
            self.app.push_screen(ConfirmScreen(f"Template save error: {e}"))
            return
        if path.exists():
            self.app.push_screen(
                ConfirmScreen(f"Template '{name}' already exists. Overwrite?"),
                lambda ok: self._do_save_template(name, overwrite=True) if ok else None,
            )
        else:
            self._do_save_template(name)

    def _do_save_template(self, name: str, overwrite: bool = False):
        state = self.app.state  # type: ignore
        try:
            path = save_user_template(name, export_protocol(state.protocol), overwrite=overwrite)
            self.notify(f"Template saved: {path}")
        except Exception as e:
            self.app.push_screen(ConfirmScreen(f"Template save error: {e}"))
