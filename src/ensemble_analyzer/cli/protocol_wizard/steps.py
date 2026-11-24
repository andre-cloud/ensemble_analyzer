"""Step-list operations of the protocol wizard (no Textual: unit-testable).

Steps reference other steps by number (``read_orbitals``, ``read_population``,
``deltascf.inherit_from``). Inserting or removing a step renumbers the steps
after it, so these references are remapped to keep pointing at the same step.

``read_orbitals`` has a default (see :func:`default_read_orbitals`). A step whose
value is the default is flagged *auto* (private key ``_ro_auto``, never saved):
its reference is recomputed after every change, so it keeps following the
previous compatible step. A value set by hand is *manual* and only renumbered.
"""
from __future__ import annotations

import copy
from typing import Any, Callable

from ensemble_analyzer.cli.protocol_wizard.fields import ML_CALCULATORS

AUTO_KEY = "_ro_auto"

Protocol = dict[str, dict[str, Any]]


def _program(step: dict) -> str | None:
    """QM program whose orbital files a step writes, None for ML/TBLite steps."""
    calc = str(step.get("calculator") or "orca").lower()
    return None if calc in ML_CALCULATORS else calc


def default_read_orbitals(previous: list[dict], step: dict) -> str:
    """Default ``read_orbitals`` of ``step`` placed after the ``previous`` steps.

    - first step, or ML/TBLite step: none;
    - ΔSCF step: the last earlier TD-DFT (``nroots``), non-ΔSCF step of the same
      program (``protocol/deltascf.py:check_protocols`` requires a TD-DFT reference);
    - otherwise the previous step, if it ran the same QM program (orbital files
      are not interchangeable between ORCA, Gaussian and NWChem).
    """
    program = _program(step)
    if not previous or program is None:
        return ""
    if step.get("deltascf"):
        for j in range(len(previous) - 1, -1, -1):
            q = previous[j]
            if q.get("nroots") and not q.get("deltascf") and _program(q) == program:
                return str(j)
        return ""
    return str(len(previous) - 1) if _program(previous[-1]) == program else ""


def orbitals_warning(previous: list[dict], step: dict, ref: str) -> str | None:
    """Why reading orbitals of ``step`` from step ``ref`` cannot work, None if it can."""
    if not ref.isdigit() or int(ref) >= len(previous):
        return None
    source = _program(previous[int(ref)])
    if source is None:
        return f"step {ref} is an ML/TBLite step: it writes no orbitals"
    if _program(step) is None:
        return "ML/TBLite steps do not read orbitals"
    if source != _program(step):
        return f"step {ref} ran {source}: its orbital file cannot be read by {_program(step)}"
    return None


# ── list <-> dict ──────────────────────────────────────────────────

def as_list(protocol: Protocol) -> list[dict]:
    return [protocol[k] for k in sorted(protocol, key=int)]


def as_dict(steps: list[dict]) -> Protocol:
    return {str(i): s for i, s in enumerate(steps)}


# ── reference remapping ────────────────────────────────────────────

def _ref_index(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value.strip())
    return None


def _remap_refs(step: dict, mapping: Callable[[int], int | None], label: str) -> tuple[dict, list[str]]:
    """Copy of ``step`` with its step references renumbered by ``mapping``.

    References whose target disappears (mapping -> None) are dropped; the
    returned warnings name them (auto ``read_orbitals`` is recomputed, no warning).
    """
    step = copy.deepcopy(step)
    lost: list[str] = []

    def remap(container: dict, key: str, name: str, warn: bool = True):
        old = _ref_index(container.get(key))
        if old is None:
            return
        new = mapping(old)
        if new is None:
            container.pop(key)
            if warn:
                lost.append(f"step {label}: {name} pointed to the removed step {old}, now unset")
        else:
            container[key] = str(new) if isinstance(container[key], str) else new

    remap(step, "read_orbitals", "read_orbitals", warn=not step.get(AUTO_KEY))
    remap(step, "read_population", "read_population")
    if isinstance(step.get("deltascf"), dict):
        remap(step["deltascf"], "inherit_from", "deltascf.inherit_from")
    if isinstance(step.get("nmr"), dict):
        remap(step["nmr"], "opt_from", "nmr.opt_from")
        remap(step["nmr"], "coupling_from", "nmr.coupling_from")
    return step, lost


# ── operations ─────────────────────────────────────────────────────

def refresh_auto(protocol: Protocol) -> Protocol:
    """Recompute ``read_orbitals`` of the auto steps."""
    steps = as_list(protocol)
    for i, step in enumerate(steps):
        if not step.get(AUTO_KEY):
            continue
        ref = default_read_orbitals(steps[:i], step)
        if ref:
            step["read_orbitals"] = ref
        else:
            step.pop("read_orbitals", None)
    return as_dict(steps)


def insert_step(protocol: Protocol, idx: int, step: dict) -> Protocol:
    """Insert ``step`` at position ``idx`` (it may only reference earlier steps)."""
    steps = as_list(protocol)
    idx = max(0, min(idx, len(steps)))
    steps = [_remap_refs(s, lambda n: n if n < idx else n + 1, str(i))[0] for i, s in enumerate(steps)]
    steps.insert(idx, copy.deepcopy(step))
    return refresh_auto(as_dict(steps))


def replace_step(protocol: Protocol, idx: int, step: dict) -> Protocol:
    steps = as_list(protocol)
    steps[idx] = copy.deepcopy(step)
    return refresh_auto(as_dict(steps))


def remove_step(protocol: Protocol, idx: int) -> tuple[Protocol, list[str]]:
    """Remove step ``idx``; returns the new protocol and the dropped references."""
    steps = as_list(protocol)
    del steps[idx]
    out, warnings = [], []
    for i, s in enumerate(steps):
        new, lost = _remap_refs(s, lambda n: n if n < idx else (None if n == idx else n - 1), str(i))
        out.append(new)
        warnings.extend(lost)
    return refresh_auto(as_dict(out)), warnings


def infer_auto(protocol: Protocol) -> Protocol:
    """Flag the steps of a loaded protocol: auto iff ``read_orbitals`` is the default.

    Steps without ``read_orbitals`` stay manual, so loading never changes a file.
    """
    steps = copy.deepcopy(as_list(protocol))
    for i, step in enumerate(steps):
        ref = str(step.get("read_orbitals") or "").strip()
        step[AUTO_KEY] = bool(ref) and ref == default_read_orbitals(steps[:i], step)
    return as_dict(steps)


def apply_defaults(protocol: Protocol) -> Protocol:
    """New protocol from a template: steps without explicit ``read_orbitals`` become auto."""
    steps = copy.deepcopy(as_list(protocol))
    for step in steps:
        step[AUTO_KEY] = not str(step.get("read_orbitals") or "").strip()
    return refresh_auto(as_dict(steps))


def export_protocol(protocol: Protocol) -> Protocol:
    """Protocol as written to JSON: wizard-private keys removed."""
    return {
        k: {kk: copy.deepcopy(v) for kk, v in protocol[k].items() if not kk.startswith("_")}
        for k in sorted(protocol, key=int)
    }
