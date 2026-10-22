import asyncio
import json
from pathlib import Path

import pytest
from textual.app import App
from textual.widgets import Input, ListView

from ensemble_analyzer.cli.protocol_wizard.main import ProtocolApp
from ensemble_analyzer.cli.protocol_wizard.screens import (
    ConfirmScreen,
    TemplateScreen,
    _InputScreen,
)
from ensemble_analyzer.cli.protocol_wizard.templates import (
    SYSTEM_TEMPLATES,
    apply,
    load_user_templates,
    save_user_template,
    user_template_path,
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


    def test_sub_folders_recursive(self, tmp_path):
        (tmp_path / "lab" / "dft").mkdir(parents=True)
        (tmp_path / "bench").mkdir()
        _write(tmp_path / "z_root.json", [{"n": 0}])
        _write(tmp_path / "lab" / "opt.json", [{"n": 1}])
        _write(tmp_path / "lab" / "dft" / "sp.json", [{"n": 2}])
        _write(tmp_path / "bench" / "dlpno.json", [{"n": 3}])
        templates, errors = load_user_templates(tmp_path)
        assert errors == []
        # top folder first, then sub-folders alphabetically
        assert list(templates) == ["z_root", "bench/dlpno", "lab/opt", "lab/dft/sp"]

    def test_hidden_entries_skipped(self, tmp_path):
        (tmp_path / ".git").mkdir()
        _write(tmp_path / ".git" / "x.json", [{}])
        _write(tmp_path / ".hidden.json", [{}])
        assert load_user_templates(tmp_path) == ({}, [])

    def test_symlinked_folder_followed(self, tmp_path):
        shared = tmp_path / "shared"
        shared.mkdir()
        _write(shared / "group_opt.json", [{"n": 1}])
        user = tmp_path / "user"
        user.mkdir()
        (user / "group").symlink_to(shared, target_is_directory=True)
        (user / "loop").symlink_to(user, target_is_directory=True)
        templates, errors = load_user_templates(user)
        assert errors == []
        assert list(templates) == ["group/group_opt"]

    def test_error_reports_relative_path(self, tmp_path):
        (tmp_path / "lab").mkdir()
        _write(tmp_path / "lab" / "broken.json", "{")
        _, errors = load_user_templates(tmp_path)
        assert errors[0].startswith("lab/broken.json:")


class TestSaveUserTemplate:

    @pytest.mark.parametrize("name, expected", [
        ("my_opt", "my_opt.json"),
        ("my_opt.json", "my_opt.json"),
        (" lab/dft/sp ", "lab/dft/sp.json"),
        ("b3lyp.d3", "b3lyp.d3.json"),
    ])
    def test_path(self, tmp_path, name, expected):
        assert user_template_path(name, tmp_path) == tmp_path / expected

    @pytest.mark.parametrize("name", ["", "  ", "/etc/x", "../x", "lab/../../x", ".hidden", "lab/.x/y"])
    def test_invalid_names(self, tmp_path, name):
        with pytest.raises(ValueError):
            user_template_path(name, tmp_path)

    def test_save_creates_folders_and_round_trips(self, tmp_path):
        protocol = {"0": {"functional": "A"}, "1": {"functional": "B"}}
        path = save_user_template("lab/new/opt", protocol, tmp_path)
        assert path == tmp_path / "lab" / "new" / "opt.json"
        templates, _ = load_user_templates(tmp_path)
        assert apply(templates["lab/new/opt"]) == protocol

    def test_overwrite(self, tmp_path):
        save_user_template("opt", {"0": {"n": 1}}, tmp_path)
        with pytest.raises(FileExistsError):
            save_user_template("opt", {"0": {"n": 2}}, tmp_path)
        save_user_template("opt", {"0": {"n": 2}}, tmp_path, overwrite=True)
        assert json.loads((tmp_path / "opt.json").read_text()) == {"0": {"n": 2}}


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

    def test_sub_folder_headers(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))
        (tmp_path / "lab").mkdir()
        _write(tmp_path / "mine.json", [{"n": 0}])
        _write(tmp_path / "lab" / "a.json", [{"n": 1}])
        _write(tmp_path / "lab" / "b.json", [{"n": 2}])

        async def check(app, pilot):
            screen = app.screen
            lv = screen.query_one("#template_list", ListView)
            labels = [str(item.query_one("Static").render()) for item in lv.children]
            user = labels[len(SYSTEM_TEMPLATES) + 1:]
            assert user[0].startswith("User Templates")
            assert user[1:] == ["  mine", "  lab/", "    a", "    b"]
            entries = screen._entries[len(SYSTEM_TEMPLATES) + 2:]
            assert entries == [[{"n": 0}], None, [{"n": 1}], [{"n": 2}]]

        self._run(check)

    def test_no_user_templates_note(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path / "missing"))

        async def check(app, pilot):
            entries = app.screen._entries
            assert len(entries) == len(SYSTEM_TEMPLATES) + 3
            assert entries[-1] is None

        self._run(check)


class TestSaveAsTemplateButton:

    def test_save_and_overwrite(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))

        async def save(pilot, name):
            await pilot.click("#save_template")
            await pilot.pause()
            assert isinstance(pilot.app.screen, _InputScreen)
            pilot.app.screen.query_one("#input_field", Input).value = name
            await pilot.click("#ok")
            await pilot.pause()

        async def go():
            app = ProtocolApp()
            async with app.run_test(size=(120, 60)) as pilot:
                await pilot.pause()
                app.state.protocol = {"0": {"functional": "A"}}
                await save(pilot, "lab/opt")
                target = tmp_path / "lab" / "opt.json"
                assert json.loads(target.read_text()) == {"0": {"functional": "A"}}

                app.state.protocol = {"0": {"functional": "B"}}
                await save(pilot, "lab/opt")
                assert isinstance(app.screen, ConfirmScreen)
                await pilot.click("#yes")
                await pilot.pause()
                assert json.loads(target.read_text()) == {"0": {"functional": "B"}}

        asyncio.run(go())

    def test_empty_protocol_not_saved(self, monkeypatch, tmp_path):
        monkeypatch.setenv("ENAN_TEMPLATES_DIR", str(tmp_path))

        async def go():
            app = ProtocolApp()
            async with app.run_test(size=(120, 60)) as pilot:
                await pilot.pause()
                await pilot.click("#save_template")
                await pilot.pause()
                assert isinstance(app.screen, ConfirmScreen)
                assert list(tmp_path.iterdir()) == []

        asyncio.run(go())
