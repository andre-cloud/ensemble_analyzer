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
        mock_protocol.functional = None
        return mock_conformer, mock_protocol

    def test_import_error_when_fairchem_missing(self, setup):
        conf, proto = setup
        with patch(
            "ensemble_analyzer.calculators.uma.FAIRChemCalculator", None
        ):
            from ensemble_analyzer.calculators.uma import UMAMlCalc
            calc = UMAMlCalc(proto, 4, conf)
            with pytest.raises(ImportError, match="fairchem-core"):
                calc._get_ml_calculator()

    @patch("ensemble_analyzer.calculators.uma.get_models_dir")
    def test_file_not_found_when_model_missing(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/uma")

        from ensemble_analyzer.calculators.uma import UMAMlCalc
        calc = UMAMlCalc(proto, 4, conf)

        with patch(
            "ensemble_analyzer.calculators.uma.FAIRChemCalculator",
            MagicMock(),
        ):
            with patch(
                "ensemble_analyzer.calculators.uma.load_predict_unit",
                MagicMock(),
            ):
                with patch(
                    "ensemble_analyzer.calculators.uma.torch"
                ) as mock_torch:
                    mock_torch.cuda.is_available.return_value = False
                    with pytest.raises(FileNotFoundError, match="Model file not found"):
                        calc._get_ml_calculator()

    @patch("ensemble_analyzer.calculators.uma.get_models_dir")
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

        with patch(
            "ensemble_analyzer.calculators.uma.FAIRChemCalculator",
            mock_fairchem,
        ):
            with patch(
                "ensemble_analyzer.calculators.uma.load_predict_unit",
                mock_load,
            ):
                with patch(
                    "ensemble_analyzer.calculators.uma.torch"
                ) as mock_torch:
                    mock_torch.cuda.is_available.return_value = False

                    result = calc._get_ml_calculator()

                    mock_load.assert_called_once_with(
                        path=model_file,
                        device="cpu",
                        inference_settings="turbo",
                    )
                    assert result is not None

    @patch("ensemble_analyzer.calculators.uma.get_models_dir")
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

        with patch(
            "ensemble_analyzer.calculators.uma.FAIRChemCalculator",
            MagicMock(),
        ):
            with patch(
                "ensemble_analyzer.calculators.uma.load_predict_unit",
                mock_load,
            ):
                with patch(
                    "ensemble_analyzer.calculators.uma.torch"
                ) as mock_torch:
                    mock_torch.cuda.is_available.return_value = False
                    calc._get_ml_calculator(method="custom_model.pt")

                    mock_load.assert_called_once_with(
                        path=custom_model,
                        device="cpu",
                        inference_settings="turbo",
                    )
