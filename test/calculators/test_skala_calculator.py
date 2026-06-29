import pytest
from unittest.mock import MagicMock, patch


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
    @pytest.fixture
    def setup(self, mock_conformer, mock_protocol):
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.functional = "skala-1.1"
        mock_protocol.basis = "def2-tzvp"
        return mock_conformer, mock_protocol

    def test_import_error_when_skala_missing(self, setup):
        conf, proto = setup
        with patch("enan_calculators._skala._Skala", None):
            from ensemble_analyzer.calculators.skala import SkalaCalc
            calc = SkalaCalc(proto, 4, conf)
            with pytest.raises(ImportError, match="skala"):
                calc._get_ml_calculator()

    def test_successful_calculator_creation(self, setup):
        conf, proto = setup
        from ensemble_analyzer.calculators.skala import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            result = calc._get_ml_calculator()
            mock_skala.assert_called_once_with(
                xc="skala-1.1", basis="def2-tzvp",
                charge=0, multiplicity=1,
            )
            assert result is mock_skala.return_value

    def test_factory_propagates_charge_and_mult(self, setup):
        conf, proto = setup
        proto.charge = -1
        proto.mult = 3

        from ensemble_analyzer.calculators.skala import SkalaCalc
        calc = SkalaCalc(proto, 4, conf)

        mock_skala = MagicMock()
        with patch("enan_calculators._skala._Skala", mock_skala):
            result = calc._get_ml_calculator()
            mock_skala.assert_called_once_with(
                xc="skala-1.1", basis="def2-tzvp",
                charge=-1, multiplicity=3,
            )
