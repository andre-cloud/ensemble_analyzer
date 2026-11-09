import asyncio

import pytest
from textual.app import App
from textual.widgets import Input

from ensemble_analyzer.cli.protocol_wizard.screens import (
    StepEditorScreen,
    _fmt_cluster,
    _parse_cluster,
    _parse_cpu_per_job,
    _step_label,
)
from ensemble_analyzer.protocol.protocol import Protocol


BASE_STEP = {"calculator": "orca", "functional": "B3LYP", "basis": "def2-SVP"}


class TestParseCpuPerJob:

    @pytest.mark.parametrize("raw", ["", "  ", "abc", "4.5", "0", "-2"])
    def test_invalid_is_unset(self, raw):
        assert _parse_cpu_per_job(raw) is None

    @pytest.mark.parametrize("raw, expected", [("8", 8), (" 1 ", 1), ("16", 16)])
    def test_valid(self, raw, expected):
        assert _parse_cpu_per_job(raw) == expected


class TestParseCluster:

    @pytest.mark.parametrize("raw", ["true", "True", "auto", "AUTO", "yes"])
    def test_true(self, raw):
        assert _parse_cluster(raw) is True

    def test_explicit_k(self):
        assert _parse_cluster("5") == 5

    @pytest.mark.parametrize("raw", ["", "false", "False", "no", "0", "1", "-3", "x"])
    def test_false(self, raw):
        assert _parse_cluster(raw) is False


class TestFmtCluster:

    @pytest.mark.parametrize("value, text", [
        (True, "true"), (False, "false"), (None, "false"), (5, "5"), (1, "false"),
    ])
    def test_fmt(self, value, text):
        assert _fmt_cluster(value) == text

    @pytest.mark.parametrize("value", [True, False, 7])
    def test_round_trip(self, value):
        assert _parse_cluster(_fmt_cluster(value)) == value


def test_step_label_shows_cpu_per_job():
    assert "[8 CPU/job]" in _step_label("0", {**BASE_STEP, "cpu_per_job": 8})
    assert "CPU/job" not in _step_label("0", BASE_STEP)


class TestStepEditorRoundTrip:

    def _read(self, step, edits=None):
        """Open the editor on ``step``, apply ``edits`` (id -> text), return _read()."""
        out = {}

        class _App(App):
            def on_mount(self):
                self.push_screen(StepEditorScreen(step=step, step_num=0, is_new=False))

        async def go():
            app = _App()
            async with app.run_test(size=(120, 60)) as pilot:
                await pilot.pause()
                screen = app.screen
                for wid, text in (edits or {}).items():
                    screen.query_one(f"#{wid}", Input).value = text
                out["step"] = screen._read()

        asyncio.run(go())
        return out["step"]

    def test_values_preserved(self):
        step = self._read({**BASE_STEP, "cpu_per_job": 8, "cluster": True})
        assert step["cpu_per_job"] == 8
        assert step["cluster"] is True
        p = Protocol(number=0, **step)
        assert p.cpu_per_job == 8
        assert p.clustering

    def test_set_from_input(self):
        step = self._read(BASE_STEP, {"f_cpu_per_job": "4", "f_cluster": "12"})
        assert step["cpu_per_job"] == 4
        assert step["cluster"] == 12
        Protocol(number=0, **step)

    @pytest.mark.parametrize("text", ["false", "", "1"])
    def test_cluster_false_not_written(self, text):
        # False is the protocol default, so _clean leaves it out of the JSON
        step = self._read({**BASE_STEP, "cluster": 5}, {"f_cluster": text})
        assert "cluster" not in step
        assert not Protocol(number=0, **step).clustering

    def test_cluster_true_from_input(self):
        step = self._read(BASE_STEP, {"f_cluster": "true"})
        assert step["cluster"] is True

    def test_blank_cpu_per_job_not_written(self):
        step = self._read({**BASE_STEP, "cpu_per_job": 8}, {"f_cpu_per_job": ""})
        assert "cpu_per_job" not in step
        Protocol(number=0, **step)
