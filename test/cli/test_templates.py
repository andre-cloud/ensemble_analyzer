import asyncio
import json
from pathlib import Path

from textual.app import App
from textual.widgets import ListView

from ensemble_analyzer.cli.protocol_wizard.screens import TemplateScreen
from ensemble_analyzer.cli.protocol_wizard.templates import (
    SYSTEM_TEMPLATES,
    apply,
    load_user_templates,
    user_templates_dir,
)


def _write(path: Path, data) -> None:
    path.write_text(json.dumps(data) if not isinstance(data, str) else data, encoding="utf-8")


class TestUserTemplatesDir:

    def test_default_dir(self, monkeypatch, tmp_path):
        monkeypatch.delenv("ENAN_TEMPLATES_DIR", raising=False)
        monkeypatch.setattr(Path, "home", lambda: tmp_path)
        assert user_templates_dir() == tmp_path / ".ensemble_analyzer" / "template"

    def test_env_override(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))
        assert user_templates_dir() == tmp_path


class TestLoadUserTemplates:

    def test_missing_dir(self, tmp_path):
        assert load_user_templates(tmp_path / "nope") == ({}, [])

    def test_protocol_format_sorted(self, tmp_path):
        _write(tmp_path / "my_opt.json", {"1": {"functional": "B"}, "0": {"functional": "A"}, "10": {"functional": "C"}})
        templates, errors = load_user_templates(tmp_path)
        assert errors == []
        assert [s["functional"] for s in templates["my_opt"]] == ["A", "B", "C"]

    def test_list_format(self, tmp_path):
        _write(tmp_path / "sp.json", [{"calculator": "orca"}])
        templates, _ = load_user_templates(tmp_path)
        assert templates == {"sp": [{"calculator": "orca"}]}

    def test_uses_env_dir_by_default(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))
        _write(tmp_path / "sp.json", [{"calculator": "orca"}])
        templates, _ = load_user_templates()
        assert list(templates) == ["sp"]

    def test_invalid_files_reported(self, tmp_path):
        _write(tmp_path / "good.json", [{"calculator": "orca"}])
        _write(tmp_path / "broken.json", "{not json")
        _write(tmp_path / "string.json", "just a string")
        _write(tmp_path / "bad_step.json", [1, 2])
        _write(tmp_path / "bad_key.json", {"a": {}})
        _write(tmp_path / "ignored.txt", "not a template")
        templates, errors = load_user_templates(tmp_path)
        assert list(templates) == ["good"]
        assert sorted(e.split(":")[0] for e in errors) == [
            "bad_key.json", "bad_step.json", "broken.json", "string.json",
        ]


class TestApply:

    def test_keys_and_deep_copy(self):
        steps = [{"constraints": [[1, 2]]}, {"opt": True}]
        protocol = apply(steps)
        assert list(protocol) == ["0", "1"]
        protocol["0"]["constraints"].append([3, 4])
        assert steps[0]["constraints"] == [[1, 2]]


class TestTemplateScreen:

    def _run(self, check):
        class _App(App):
            def on_mount(self):
                self.result = "unset"
                self.push_screen(TemplateScreen(), lambda r: setattr(self, "result", r))

        async def go():
            app = _App()
            async with app.run_test() as pilot:
                await pilot.pause()
                await check(app, pilot)

        asyncio.run(go())

    def test_sections_and_selection(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))
        _write(tmp_path / "mine.json", {"0": {"functional": "PBE0"}})
        _write(tmp_path / "broken.json", "{")

        async def check(app, pilot):
            screen = app.screen
            lv = screen.query_one("#template_list", ListView)
            entries = screen._entries
            n_sys = len(SYSTEM_TEMPLATES)
            # header, system templates, header, user template, error note
            assert len(entries) == n_sys + 4
            assert entries[0] is None and entries[n_sys + 1] is None
            assert entries[n_sys + 2] == [{"functional": "PBE0"}]
            assert entries[-1] is None
            assert lv.children[0].disabled and lv.children[-1].disabled
            lv.index = n_sys + 2
            await pilot.press("enter")
            await pilot.pause()
            assert app.result == [{"functional": "PBE0"}]

        self._run(check)

    def test_no_user_templates_note(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path / "missing"))

        async def check(app, pilot):
            entries = app.screen._entries
            assert len(entries) == len(SYSTEM_TEMPLATES) + 3
            assert entries[-1] is None

        self._run(check)
