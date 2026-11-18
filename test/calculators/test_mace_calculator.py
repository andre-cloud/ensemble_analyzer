import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path
from ase.calculators.calculator import all_changes


class TestMACERegistration:
    def test_registered(self):
        from ensemble_analyzer.calculators import CALCULATOR_REGISTRY
        assert "mace" in CALCULATOR_REGISTRY
        assert CALCULATOR_REGISTRY["mace"].__name__ == "MACEMlCalc"

    def test_in_ml_calculators(self):
        from ensemble_analyzer.calculators.base import ML_CALCULATORS
        assert "mace" in ML_CALCULATORS

    def test_in_regex_parsing(self):
        from ensemble_analyzer.constants import regex_parsing
        assert "mace" in regex_parsing
        assert regex_parsing["mace"]["ext"] is None


class TestMACEMlCalc:
    @pytest.fixture
    def setup(self, mock_conformer, mock_protocol):
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.functional = "MACE_model.pt"
        return mock_conformer, mock_protocol

    def test_import_error_when_mace_missing(self, setup):
        conf, proto = setup
        with patch(
            "enan_calculators._mace.MACECalculator", None
        ):
            from ensemble_analyzer.calculators.mace import MACEMlCalc
            calc = MACEMlCalc(proto, 4, conf)
            with pytest.raises(ImportError, match="mace"):
                calc._get_ml_calculator()

    @patch("enan_calculators._mace.get_models_dir")
    def test_file_not_found_when_model_missing(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/mace")

        from ensemble_analyzer.calculators.mace import MACEMlCalc
        calc = MACEMlCalc(proto, 4, conf)

        with patch(
            "enan_calculators._mace.MACECalculator",
            MagicMock(),
        ):
            with patch(
                "enan_calculators._mace.torch"
            ) as mock_torch:
                mock_torch.cuda.is_available.return_value = False
                with pytest.raises(FileNotFoundError, match="MACE model not found"):
                    calc._get_ml_calculator()

    @patch("enan_calculators._mace.get_models_dir")
    def test_successful_calculator_creation(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        model_dir = tmp_path / "models" / "mace"
        model_dir.mkdir(parents=True)
        model_file = model_dir / "MACE_model.pt"
        model_file.touch()

        mock_get_models_dir.return_value = model_dir

        from ensemble_analyzer.calculators.mace import MACEMlCalc
        calc = MACEMlCalc(proto, 4, conf)

        mock_mace = MagicMock()
        mock_mace_instance = MagicMock()
        mock_mace_instance.implemented_properties = ["energy", "forces"]

        with patch(
            "enan_calculators._mace.MACECalculator",
            mock_mace,
        ):
            with patch(
                "enan_calculators._mace.torch"
            ) as mock_torch:
                mock_torch.cuda.is_available.return_value = False

                result = calc._get_ml_calculator()

                mock_mace.assert_called_once_with(
                    model_paths=str(model_file),
                    device="cpu",
                    default_dtype="float64",
                )

    @patch("enan_calculators._mace.get_models_dir")
    def test_calculator_uses_method_from_kwargs(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        model_dir = tmp_path / "models" / "mace"
        model_dir.mkdir(parents=True)
        custom_model = model_dir / "custom_model.pt"
        custom_model.touch()

        mock_get_models_dir.return_value = model_dir

        from ensemble_analyzer.calculators.mace import MACEMlCalc
        calc = MACEMlCalc(proto, 4, conf)

        mock_mace = MagicMock()
        mock_mace_instance = MagicMock()
        mock_mace_instance.implemented_properties = ["energy", "forces"]

        with patch(
            "enan_calculators._mace.MACECalculator",
            mock_mace,
        ):
            with patch(
                "enan_calculators._mace.torch"
            ) as mock_torch:
                mock_torch.cuda.is_available.return_value = False
                calc._get_ml_calculator(method="custom_model.pt")

                mock_mace.assert_called_once_with(
                    model_paths=str(custom_model),
                    device="cpu",
                    default_dtype="float64",
                )

    @patch("enan_calculators._mace.get_models_dir")
    def test_foundation_model(self, mock_get_models_dir, setup):
        conf, proto = setup

        mock_result = MagicMock(implemented_properties=["energy", "forces"])

        with patch(
            "enan_calculators._mace.create_mace_calc",
            return_value=mock_result,
        ):
            from ensemble_analyzer.calculators.mace import MACEMlCalc
            calc = MACEMlCalc(proto, 4, conf)
            result = calc._get_ml_calculator(method="mp")
            assert result is not None

    def test_foundation_mapping(self):
        from enan_calculators._mace import _FOUNDATION
        assert _FOUNDATION["mp"] == "mace_mp"
        assert _FOUNDATION["off"] == "mace_off"
        assert _FOUNDATION["anicc"] == "mace_anicc"
        assert _FOUNDATION["mdp"] == "mace_mdp"
        assert _FOUNDATION["omol"] == "mace_omol"

    def test_missing_method_raises(self, setup):
        conf, proto = setup
        proto.functional = None
        from ensemble_analyzer.calculators.mace import MACEMlCalc
        calc = MACEMlCalc(proto, 4, conf)
        with pytest.raises(ValueError, match="requires a model"):
            calc._get_ml_calculator()

    @pytest.mark.parametrize("method, charge, mult, warns", [
        ("off", 1, 1, True), ("mp", 0, 2, True), ("off", 0, 1, False), ("omol", 1, 2, False),
    ])
    def test_charge_spin_warning(self, method, charge, mult, warns):
        import warnings
        import enan_calculators._mace as m
        with patch.object(m, "MACECalculator", MagicMock()), \
             patch.object(m, "torch") as mock_torch, \
             patch.object(m, "_load_foundation", return_value=MagicMock(implemented_properties=["energy"])), \
             patch.object(m, "_NEEDS_REBUILD", False), \
             patch.dict(m._PREDICTOR_CACHE, clear=True):
            mock_torch.cuda.is_available.return_value = False
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                m.create_mace_calc(charge, mult, method)
        assert any("ignores charge and spin" in str(w.message) for w in rec) is warns

    def test_wrapped_calc_injects_charge_and_spin(self, setup):
        conf, proto = setup
        proto.charge = -1
        proto.mult = 3

        from enan_calculators._mace import _MACEWrappedCalc

        mock_inner = MagicMock()
        mock_inner.implemented_properties = ["energy", "forces"]
        mock_inner.results = {"energy": -1.0}

        wrapped = _MACEWrappedCalc(mock_inner, proto.charge, proto.mult)

        mock_atoms = MagicMock()
        mock_atoms.info = {}

        wrapped.calculate(mock_atoms)
        assert mock_atoms.info["charge"] == -1
        assert mock_atoms.info["spin"] == 3
        mock_inner.calculate.assert_called_once_with(mock_atoms, None, all_changes)
