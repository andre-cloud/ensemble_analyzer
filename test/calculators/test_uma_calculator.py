import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path


class TestUMARegistration:
    def test_registered(self):
        from ensemble_analyzer.calculators import CALCULATOR_REGISTRY
        assert "uma" in CALCULATOR_REGISTRY
        assert CALCULATOR_REGISTRY["uma"].__name__ == "UMAMlCalc"

    def test_in_ml_calculators(self):
        from ensemble_analyzer.calculators.base import ML_CALCULATORS
        assert "uma" in ML_CALCULATORS


class TestUMAMlCalc:
    @pytest.fixture
    def setup(self, mock_conformer, mock_protocol):
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.functional = "uma-s-1.pt"
        return mock_conformer, mock_protocol

    def test_requires_method(self, setup):
        conf, proto = setup
        proto.functional = None
        from ensemble_analyzer.calculators.uma import UMAMlCalc
        calc = UMAMlCalc(proto, 4, conf)
        with pytest.raises(ValueError, match="requires a model path"):
            calc._get_ml_calculator()

    def test_import_error_when_fairchem_missing(self, setup):
        conf, proto = setup
        with patch("enan_calculators._ml_inference.FAIRChemCalculator", None):
            from ensemble_analyzer.calculators.uma import UMAMlCalc
            calc = UMAMlCalc(proto, 4, conf)
            with pytest.raises(ImportError, match="fairchem-core"):
                calc._get_ml_calculator()

    @patch("enan_calculators._ml_inference.get_models_dir")
    def test_file_not_found_when_model_missing(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/uma")

        from ensemble_analyzer.calculators.uma import UMAMlCalc
        calc = UMAMlCalc(proto, 4, conf)

        with patch("enan_calculators._ml_inference.FAIRChemCalculator", MagicMock()):
            with patch("enan_calculators._ml_inference.load_predict_unit", MagicMock()):
                with patch("enan_calculators._ml_inference.torch") as mock_torch:
                    mock_torch.cuda.is_available.return_value = False
                    with pytest.raises(FileNotFoundError, match="Model file not found"):
                        calc._get_ml_calculator()

    @patch("enan_calculators._ml_inference.get_models_dir")
    def test_successful_calculator_creation(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        model_dir = tmp_path / "models" / "uma"
        model_dir.mkdir(parents=True)
        model_file = model_dir / "uma-s-1.pt"
        model_file.touch()

        mock_get_models_dir.return_value = model_dir

        from ensemble_analyzer.calculators.uma import UMAMlCalc
        calc = UMAMlCalc(proto, 4, conf)

        mock_fairchem = MagicMock()
        mock_predictor = MagicMock()
        mock_load = MagicMock(return_value=mock_predictor)

        with patch("enan_calculators._ml_inference.FAIRChemCalculator", mock_fairchem):
            with patch("enan_calculators._ml_inference.load_predict_unit", mock_load):
                with patch("enan_calculators._ml_inference.torch") as mock_torch:
                    mock_torch.cuda.is_available.return_value = False

                    result = calc._get_ml_calculator()

                    mock_load.assert_called_once()
                    assert mock_load.call_args[1]["path"] == model_file
                    assert mock_load.call_args[1]["device"] == "cpu"
                    assert "inference_settings" in mock_load.call_args[1]
                    assert result is not None

    @patch("enan_calculators._ml_inference.get_models_dir")
    def test_calculator_uses_method_from_kwargs(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        model_dir = tmp_path / "models" / "uma"
        model_dir.mkdir(parents=True)
        custom_model = model_dir / "custom_model.pt"
        custom_model.touch()

        mock_get_models_dir.return_value = model_dir

        from ensemble_analyzer.calculators.uma import UMAMlCalc
        calc = UMAMlCalc(proto, 4, conf)

        mock_predictor = MagicMock()
        mock_load = MagicMock(return_value=mock_predictor)

        with patch("enan_calculators._ml_inference.FAIRChemCalculator", MagicMock()):
            with patch("enan_calculators._ml_inference.load_predict_unit", mock_load):
                with patch("enan_calculators._ml_inference.torch") as mock_torch:
                    mock_torch.cuda.is_available.return_value = False
                    calc._get_ml_calculator(method="custom_model.pt")

                    mock_load.assert_called_once()
                    assert mock_load.call_args[1]["path"] == custom_model
                    assert mock_load.call_args[1]["device"] == "cpu"
                    assert "inference_settings" in mock_load.call_args[1]
