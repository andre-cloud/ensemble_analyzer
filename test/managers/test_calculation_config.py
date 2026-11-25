import pytest
import json
from unittest.mock import MagicMock, patch, mock_open
from pathlib import Path
from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer.protocol.protocol import Protocol


def _args(**over):
    from types import SimpleNamespace
    d = dict(cpu=1, temperature=298.15, definition=4, fwhm_vibro=None,
             fwhm_electro=None, shift_vibro=None, shift_electro=None,
             interest_vibro=None, interest_electro=None, invert=False,
             include_H=False, restart=False, pressure=101.325, linear=False,
             alpha=4, cut_off=100, norm_exp=False)
    d.update(over)
    return SimpleNamespace(**d)


class TestCalculationConfig:

    def test_default_values(self):
        cfg = CalculationConfig()
        assert cfg.cpu == 1
        assert cfg.temperature == 298.15
        assert cfg.definition == 4
        assert cfg.fwhm == {"vibro": None, "electro": None}

    def test_post_init_defaults(self):
        cfg = CalculationConfig(cpu=8)
        assert cfg.fwhm == {"vibro": None, "electro": None}
        assert cfg.shift == {"vibro": None, "electro": None}
        assert cfg.interested == {"vibro": None, "electro": None}

    def test_from_args_no_existing_settings(self):
        args = _args(cpu=4, temperature=300.0, definition=3, restart=True)

        m = mock_open()
        with patch("ensemble_analyzer._managers.calculation_config.Path.exists", return_value=False):
            with patch("builtins.open", m):
                cfg = CalculationConfig.from_args(args)
                assert cfg.cpu == 4

    def test_from_args_with_existing_settings(self):
        args = _args()

        settings_data = json.dumps({"cpu": 16, "temperature": 350.0})
        m = mock_open(read_data=settings_data)
        with patch("ensemble_analyzer._managers.calculation_config.Path.exists", return_value=True):
            with patch("builtins.open", m):
                cfg = CalculationConfig.from_args(args)
                assert cfg.cpu == 16
                assert cfg.temperature == 350.0

    def test_load(self, tmp_path):
        fp = tmp_path / "settings.json"
        fp.write_text(json.dumps({"cpu": 8, "temperature": 300.0, "definition": 3, "fwhm_vibro": 5}))
        loaded = CalculationConfig.load(fp)
        assert loaded.cpu == 8
        assert loaded.temperature == 300.0
        assert loaded.definition == 3
        assert loaded.fwhm == {"vibro": 5, "electro": None}

    def test_from_args_thermo_flags(self, tmp_path, monkeypatch):
        """--pressure, --linear, --alpha and --cut-off reach the config."""
        monkeypatch.chdir(tmp_path)
        args = _args(pressure=50.0, linear=True, alpha=3.0, cut_off=80.0)
        cfg = CalculationConfig.from_args(args)
        assert (cfg.P, cfg.linear, cfg.alpha, cfg.cut_off) == (50.0, True, 3.0, 80.0)
        assert json.loads((tmp_path / "settings.json").read_text())["P"] == 50.0

    def test_load_missing(self):
        with pytest.raises(FileNotFoundError):
            CalculationConfig.load(Path("/nonexistent/settings.json"))

    def test_create_log(self, mock_protocol):
        cfg = CalculationConfig(cpu=4, temperature=298.15)
        result = cfg.create_log([mock_protocol], 10)
        assert result["conformers"] == 10
        assert result["len_protocols"] == 1
        assert result["cpu"] == 4
