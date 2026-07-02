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

    def test_import_error_when_skala_missing(self, setup):
        conf, proto = setup
        with patch("enan_calculators._skala._load_skala",
                   side_effect=ImportError("skala module missing")):
            from ensemble_analyzer.calculators.skala_calc import SkalaCalc
            calc = SkalaCalc(proto, 4, conf)
            with pytest.raises(ImportError, match="skala"):
                calc._get_ml_calculator()

    @patch("enan_calculators._skala.get_models_dir")
    def test_local_model_discovered(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        models_dir = tmp_path / "models" / "skala"
        models_dir.mkdir(parents=True)
        model_file = models_dir / "skala-1.1"
        model_file.touch()
        mock_get_models_dir.return_value = models_dir

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            result = calc._get_ml_calculator()
            assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(model_file)
            mock_skala.assert_called_once_with(
                xc="skala-1.1", basis="def2-tzvp",
                charge=0, multiplicity=1,
            )

    @patch("enan_calculators._skala.get_models_dir")
    def test_local_model_with_fun_extension(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        models_dir = tmp_path / "models" / "skala"
        models_dir.mkdir(parents=True)
        model_file = models_dir / "skala-1.1.fun"
        model_file.touch()
        mock_get_models_dir.return_value = models_dir

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            calc._get_ml_calculator()
            assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(model_file)

    @patch("enan_calculators._skala.get_models_dir")
    def test_no_local_model_skips_env_var(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/skala")

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            assert "SKALA_LOCAL_MODEL_PATH" not in os.environ
            calc._get_ml_calculator()

    @patch("enan_calculators._skala.get_models_dir")
    def test_clears_stale_env_var(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/skala")

        os.environ["SKALA_LOCAL_MODEL_PATH"] = "/stale/path"

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            calc._get_ml_calculator()
            assert "SKALA_LOCAL_MODEL_PATH" not in os.environ

    @patch("enan_calculators._skala.get_models_dir")
    def test_direct_path_fallback(self, mock_get_models_dir, setup, tmp_path):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/skala")

        custom_path = tmp_path / "custom_model.fun"
        custom_path.touch()
        proto.functional = str(custom_path)

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            calc._get_ml_calculator()
            assert os.environ["SKALA_LOCAL_MODEL_PATH"] == str(custom_path)

    @patch("enan_calculators._skala.get_models_dir")
    def test_cache_reuses_same_config(self, mock_get_models_dir, setup):
        conf, proto = setup
        mock_get_models_dir.return_value = Path("/nonexistent/models/skala")

        from ensemble_analyzer.calculators.skala_calc import SkalaCalc
        c1 = SkalaCalc(proto, 4, conf)
        c2 = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            assert c1._get_ml_calculator() is c2._get_ml_calculator()
            mock_skala.assert_called_once()
