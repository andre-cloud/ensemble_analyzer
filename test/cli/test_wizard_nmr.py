import pytest

from ensemble_analyzer.cli.protocol_wizard.nmr_fields import (
    format_decouple,
    format_groups,
    format_references,
    parse_decouple,
    parse_groups,
    parse_indices,
    parse_references,
)
from ensemble_analyzer.cli.protocol_wizard.steps import _remap_refs


def test_references_round_trip():
    text = "TMS: 1H 13C; PhF = Fc1ccccc1: 13C; mine = refs/cfcl3.xyz: 19F"
    refs = parse_references(text)
    assert refs == [
        {"name": "TMS", "nuclei": ["1H", "13C"]},
        {"name": "PhF", "smiles": "Fc1ccccc1", "nuclei": ["13C"]},
        {"name": "mine", "xyz": "refs/cfcl3.xyz", "nuclei": ["19F"]},
    ]
    assert parse_references(format_references(refs)) == refs


@pytest.mark.parametrize("bad", ["TMS", "TMS:", ": 1H"])
def test_references_errors(bad):
    with pytest.raises(ValueError):
        parse_references(bad)


def test_decouple_groups_and_indices():
    assert parse_decouple("13C: 1H; 15N: 1H 13C") == {"13C": ["1H"], "15N": ["1H", "13C"]}
    assert parse_decouple(format_decouple({"13C": ["1H"]})) == {"13C": ["1H"]}
    assert parse_indices("0, 3 4") == [0, 3, 4] and parse_indices("") is None
    assert parse_groups(format_groups([[1, 2, 3], [7, 8]])) == [[1, 2, 3], [7, 8]]


def test_remap_nmr_step_references():
    step = {"nmr": {"opt_from": 2, "coupling_from": 3}}
    out, lost = _remap_refs(step, lambda i: None if i == 3 else i - 1, "4")
    assert out["nmr"]["opt_from"] == 1 and "coupling_from" not in out["nmr"] and lost
    out, _ = _remap_refs({"nmr": {"coupling_from": False}}, lambda i: i + 1, "1")
    assert out["nmr"]["coupling_from"] is False


class TestEditorNMR:
    """The NMR section round-trips through the editor and gives a valid Protocol."""

    def _read(self, step, edits=None):
        import asyncio
        from textual.app import App
        from textual.widgets import Input
        from ensemble_analyzer.cli.protocol_wizard.screens import StepEditorScreen

        out = {}

        class _App(App):
            def on_mount(self):
                self.push_screen(StepEditorScreen(step=step, step_num=1, is_new=False))

        async def go():
            app = _App()
            async with app.run_test(size=(120, 60)) as pilot:
                await pilot.pause()
                for wid, text in (edits or {}).items():
                    app.screen.query_one(f"#{wid}", Input).value = text
                out["step"] = {k: v for k, v in app.screen._read().items() if not k.startswith("_")}

        asyncio.run(go())
        return out["step"]

    def test_round_trip(self):
        from ensemble_analyzer.protocol.protocol import Protocol
        nmr = {"shielding": True, "coupling": True, "nuclei": ["1H", "13C"], "atoms": [0, 3],
               "references": [{"name": "TMS", "nuclei": ["1H", "13C"]}], "decouple": {"13C": ["1H"]},
               "equivalent": [[4, 5]], "coupling_from": False, "opt_from": 0}
        step = self._read({"calculator": "orca", "functional": "PBE0", "basis": "pcSseg-1", "nmr": nmr})
        p = Protocol(number=1, **step)
        assert p.nmr["nuclei"] == ["1H", "13C"] and p.nmr["atoms"] == [0, 3]
        assert p.nmr["decouple"] == {"13C": ["1H"]} and p.nmr["equivalent"] == [[4, 5]]
        assert p.nmr["coupling_from"] is False and p.nmr["opt_from"] == 0

    def test_disabled_step_has_no_nmr(self):
        step = self._read({"calculator": "orca", "functional": "PBE0", "basis": "def2-SVP"})
        assert "nmr" not in step
