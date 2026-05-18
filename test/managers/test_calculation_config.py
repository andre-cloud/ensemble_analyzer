import pytest
import json
from unittest.mock import MagicMock, patch, mock_open
from pathlib import Path
from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._protocol.protocol import Protocol


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

    def test_validate_ok(self):
        cfg = CalculationConfig(cpu=4, temperature=300.0, definition=3)
        cfg.validate()

    def test_validate_bad_cpu(self):
        cfg = CalculationConfig(cpu=0)
        with pytest.raises(ValueError, match="CPU count must be ≥ 1"):
            cfg.validate()

    def test_validate_bad_temperature(self):
        cfg = CalculationConfig(temperature=-1)
        with pytest.raises(ValueError, match="Temperature must be > 0"):
            cfg.validate()

    def test_validate_bad_definition(self):
        cfg = CalculationConfig(definition=0)
        with pytest.raises(ValueError, match="Definition must be ≥ 1"):
            cfg.validate()

    def test_from_args_no_existing_settings(self):
        args = MagicMock()
        args.cpu = 4
        args.temperature = 300.0
        args.definition = 3
        args.fwhm_vibro = None
        args.fwhm_electro = None
        args.shift_vibro = None
        args.shift_electro = None
        args.interest_vibro = None
        args.interest_electro = None
        args.invert = False
        args.exclude_H = False
        args.restart = True

        m = mock_open()
        with patch("ensemble_analyzer._managers.calculation_config.Path.exists", return_value=False):
            with patch("builtins.open", m):
                cfg = CalculationConfig.from_args(args)
                assert cfg.cpu == 4

    def test_from_args_with_existing_settings(self):
        args = MagicMock()
        args.cpu = 1
        args.temperature = 298.15
        args.definition = 4
        args.fwhm_vibro = None
        args.fwhm_electro = None
        args.shift_vibro = None
        args.shift_electro = None
        args.interest_vibro = None
        args.interest_electro = None
        args.invert = False
        args.exclude_H = False
        args.restart = False

        settings_data = json.dumps({"cpu": 16, "temperature": 350.0})
        m = mock_open(read_data=settings_data)
        with patch("ensemble_analyzer._managers.calculation_config.Path.exists", return_value=True):
            with patch("builtins.open", m):
                cfg = CalculationConfig.from_args(args)
                assert cfg.cpu == 16
                assert cfg.temperature == 350.0

    def test_to_dict(self):
        cfg = CalculationConfig(cpu=8, temperature=300.0)
        d = cfg.to_dict()
        assert d["cpu"] == 8
        assert d["temperature"] == 300.0

    def test_save_and_load(self, tmp_path):
        cfg = CalculationConfig(cpu=8, temperature=300.0, definition=3)
        fp = tmp_path / "settings.json"
        cfg.save(fp)
        assert fp.exists()

        loaded = CalculationConfig.load(fp)
        assert loaded.cpu == 8
        assert loaded.temperature == 300.0
        assert loaded.definition == 3

    def test_load_missing(self):
        with pytest.raises(FileNotFoundError):
            CalculationConfig.load(Path("/nonexistent/settings.json"))

    def test_create_log(self, mock_protocol):
        cfg = CalculationConfig(cpu=4, temperature=298.15)
        result = cfg.create_log([mock_protocol], 10)
        assert result["conformers"] == 10
        assert result["len_protocols"] == 1
        assert result["cpu"] == 4
