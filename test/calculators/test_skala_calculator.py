import os
import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path


def _clear_cache():
    from enan_calculators._skala import _PREDICTOR_CACHE
    _PREDICTOR_CACHE.clear()


class TestSkalaRegistration:
    def test_registered(self):
        from ensemble_analyzer.calculators import CALCULATOR_REGISTRY
        assert "skala" in CALCULATOR_REGISTRY
        assert CALCULATOR_REGISTRY["skala"].__name__ == "SkalaCalc"

    def test_in_ml_calculators(self):
        from ensemble_analyzer.calculators.base import ML_CALCULATORS
        assert "skala" in ML_CALCULATORS

    def test_in_regex_parsing(self):
        from ensemble_analyzer.constants import regex_parsing
        assert "skala" in regex_parsing
        assert regex_parsing["skala"]["ext"] is None


class TestSkalaCalc:
    @pytest.fixture(autouse=True)
    def clear_cache(self):
        _clear_cache()
        os.environ.pop("SKALA_LOCAL_MODEL_PATH", None)
        yield

    @pytest.fixture
    def setup(self, mock_conformer, mock_protocol):
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.functional = "skala-1.1"
        mock_protocol.basis = "def2-tzvp"
        return mock_conformer, mock_protocol

    def _calc(self, proto, conf):
        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        return SkalaCalc(proto, 4, conf)

    def test_import_error_when_skala_missing(self, setup):
        conf, proto = setup
        with patch("enan_calculators._skala.Skala", None):
            with pytest.raises(ImportError, match="skala"):
                self._calc(proto, conf)._get_ml_calculator()

    @patch("enan_calculators._skala.get_models_dir")
    def test_local_model_discovered(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = tmp_path
        model_file = tmp_path / "skala-1.1"
        model_file.touch()

        mock_skala = MagicMock()
        with patch("enan_calculators._skala.Skala", mock_skala):
            self._calc(proto, conf)._get_ml_calculator()
        assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(model_file)
        mock_skala.assert_called_once_with(
            ks_config={"functional_path": model_file},
            xc="skala-1.1", basis="def2-tzvp",
            charge=0, multiplicity=1, with_density_fit=True,
        )

    @patch("enan_calculators._skala.get_models_dir")
    def test_local_model_with_fun_extension(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = tmp_path
        model_file = tmp_path / "skala-1.1.fun"
        model_file.touch()

        mock_skala = MagicMock()
        with patch("enan_calculators._skala.Skala", mock_skala):
            self._calc(proto, conf)._get_ml_calculator()
        assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(model_file)
        assert mock_skala.call_args.kwargs["xc"] == "skala-1.1"

    @patch("enan_calculators._skala.get_models_dir")
    def test_missing_model_raises(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = tmp_path
        with patch("enan_calculators._skala.Skala", MagicMock()):
            with pytest.raises(FileNotFoundError, match="Skala model not found"):
                self._calc(proto, conf)._get_ml_calculator()

    @patch("enan_calculators._skala.get_models_dir")
    def test_direct_path_fallback(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = tmp_path / "empty"
        custom_path = tmp_path / "custom_model.fun"
        custom_path.touch()
        proto.functional = str(custom_path)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala.Skala", mock_skala):
            self._calc(proto, conf)._get_ml_calculator()
        assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(custom_path)

    @patch("enan_calculators._skala.get_models_dir")
    def test_cache_reuses_same_config(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = tmp_path
        (tmp_path / "skala-1.1").touch()

        mock_skala = MagicMock()
        with patch("enan_calculators._skala.Skala", mock_skala):
            c1 = self._calc(proto, conf)._get_ml_calculator()
            c2 = self._calc(proto, conf)._get_ml_calculator()
        assert c1 is c2
        mock_skala.assert_called_once()
