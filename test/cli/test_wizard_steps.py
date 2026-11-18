import pytest

from ensemble_analyzer.cli.protocol_wizard.steps import (
    AUTO_KEY,
    apply_defaults,
    default_read_orbitals,
    export_protocol,
    infer_auto,
    insert_step,
    orbitals_warning,
    remove_step,
    replace_step,
)

ORCA = {"calculator": "orca", "functional": "B3LYP", "basis": "def2-SVP"}
GAU = {"calculator": "gaussian", "functional": "B3LYP", "basis": "def2-SVP"}
XTB = {"calculator": "tblite", "functional": "GFN2-xTB"}
UMA = {"calculator": "uma", "functional": "uma-s-1"}
TDDFT = {**ORCA, "nroots": 5}
DSCF = {**ORCA, "deltascf": {"excitation": "auto"}}


def auto(step):
    return {**step, AUTO_KEY: True}


class TestDefault:

    def test_first_step_none(self):
        assert default_read_orbitals([], ORCA) == ""

    def test_previous_same_program(self):
        assert default_read_orbitals([ORCA, ORCA], ORCA) == "1"

    @pytest.mark.parametrize("current", [XTB, UMA])
    def test_ml_current_none(self, current):
        assert default_read_orbitals([ORCA], current) == ""

    @pytest.mark.parametrize("prev", [XTB, UMA, GAU])
    def test_incompatible_previous_none(self, prev):
        assert default_read_orbitals([ORCA, prev], ORCA) == ""

    def test_deltascf_picks_last_tddft(self):
        assert default_read_orbitals([ORCA, TDDFT, ORCA], DSCF) == "1"

    def test_deltascf_skips_deltascf_and_other_program(self):
        steps = [TDDFT, {**GAU, "nroots": 5}, {**TDDFT, "deltascf": {"excitation": "auto"}}]
        assert default_read_orbitals(steps, DSCF) == "0"

    def test_deltascf_without_tddft_none(self):
        assert default_read_orbitals([ORCA], DSCF) == ""


class TestInsert:

    def test_auto_follows_previous_step(self):
        p = apply_defaults({"0": ORCA, "1": ORCA, "2": ORCA})
        assert p["2"]["read_orbitals"] == "1"
        p = insert_step(p, 2, auto({**ORCA, "comment": "new"}))
        assert p["2"]["comment"] == "new"
        assert p["2"]["read_orbitals"] == "1"
        assert p["3"]["read_orbitals"] == "2"

    def test_manual_follows_original_step(self):
        p = {"0": ORCA, "1": ORCA, "2": {**ORCA, "read_orbitals": "1", AUTO_KEY: False}}
        p = insert_step(p, 1, auto(ORCA))
        assert p["3"]["read_orbitals"] == "2"

    def test_inserting_ml_step_breaks_the_chain(self):
        p = apply_defaults({"0": ORCA, "1": ORCA})
        p = insert_step(p, 1, auto(XTB))
        assert "read_orbitals" not in p["1"]
        assert "read_orbitals" not in p["2"]

    def test_other_refs_shift(self):
        p = {"0": ORCA, "1": TDDFT, "2": {**ORCA, "read_population": "1",
                                          "deltascf": {"excitation": "inherit", "inherit_from": 1}}}
        p = insert_step(p, 0, auto(ORCA))
        assert p["3"]["read_population"] == "2"
        assert p["3"]["deltascf"]["inherit_from"] == 2

    def test_append(self):
        p = apply_defaults({"0": ORCA})
        p = insert_step(p, 1, auto(ORCA))
        assert list(p) == ["0", "1"] and p["1"]["read_orbitals"] == "0"


class TestRemove:

    def test_shift_and_auto_recompute(self):
        p = apply_defaults({"0": ORCA, "1": ORCA, "2": ORCA, "3": ORCA})
        p, lost = remove_step(p, 1)
        assert [p[k].get("read_orbitals") for k in p] == [None, "0", "1"]
        assert lost == []

    def test_manual_ref_to_removed_step_dropped_with_warning(self):
        p = {"0": ORCA, "1": ORCA, "2": {**ORCA, "read_orbitals": "1", "read_population": "1"}}
        p, lost = remove_step(p, 1)
        assert "read_orbitals" not in p["1"] and "read_population" not in p["1"]
        assert len(lost) == 2

    def test_manual_ref_after_removed_shifts(self):
        p = {"0": ORCA, "1": ORCA, "2": ORCA, "3": {**ORCA, "read_orbitals": "2"}}
        p, _ = remove_step(p, 0)
        assert p["2"]["read_orbitals"] == "1"


class TestLoadAndExport:

    def test_infer_auto_does_not_change_steps(self):
        raw = {"0": ORCA, "1": ORCA, "2": {**ORCA, "read_orbitals": "1"}, "3": {**ORCA, "read_orbitals": "0"}}
        p = infer_auto(raw)
        assert export_protocol(p) == raw
        assert [p[k][AUTO_KEY] for k in p] == [False, False, True, False]

    def test_export_strips_private_keys(self):
        p = apply_defaults({"0": ORCA, "1": ORCA})
        out = export_protocol(p)
        assert all(AUTO_KEY not in s for s in out.values())
        assert out["1"]["read_orbitals"] == "0"

    def test_template_explicit_ref_kept(self):
        p = apply_defaults({"0": ORCA, "1": TDDFT, "2": {**DSCF, "read_orbitals": "1"}})
        assert p["2"]["read_orbitals"] == "1" and not p["2"][AUTO_KEY]

    def test_replace_step_updates_following_auto(self):
        p = apply_defaults({"0": ORCA, "1": ORCA, "2": ORCA})
        p = replace_step(p, 1, auto(XTB))
        assert "read_orbitals" not in p["1"] and "read_orbitals" not in p["2"]


def test_orbitals_warning():
    assert orbitals_warning([ORCA, XTB], ORCA, "1")
    assert orbitals_warning([GAU], ORCA, "0")
    assert orbitals_warning([ORCA], ORCA, "0") is None
