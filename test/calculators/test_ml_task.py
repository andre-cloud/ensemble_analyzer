import pytest
from unittest.mock import MagicMock, patch


class TestProtocolTask:
    def _protocol(self, **kwargs):
        from ensemble_analyzer.protocol.protocol import Protocol
        kwargs.setdefault("functional", "model.pt")
        return Protocol(number=1, **kwargs)

    def test_default_is_none(self):
        assert self._protocol(calculator="uma").task is None

    def test_value_is_kept(self):
        assert self._protocol(calculator="uma", task="omat").task == "omat"

    @pytest.mark.parametrize("value", ["", "   "])
    def test_blank_becomes_none(self, value):
        assert self._protocol(calculator="mace", task=value).task is None

    def test_aimnet_with_task_raises(self):
        with pytest.raises(ValueError, match="has no task/head"):
            self._protocol(calculator="aimnet", task="omol")

    def test_aimnet_without_task_is_fine(self):
        assert self._protocol(calculator="aimnet").task is None

    def test_qm_engine_ignores_task(self):
        assert self._protocol(calculator="orca", functional="B3LYP", task="omol").task == "omol"


class TestFAIRChemTask:
    @pytest.mark.parametrize("engine, module, factory", [
        ("uma", "enan_calculators._uma", "create_uma_calc"),
        ("fairchem", "enan_calculators._fairchem", "create_fairchem_calc"),
    ])
    @pytest.mark.parametrize("task, expected", [(None, "omol"), ("omat", "omat")])
    def test_task_name(self, engine, module, factory, task, expected):
        import importlib
        mod = importlib.import_module(module)
        with patch(f"{module}.create_ml_calc") as mock_create:
            getattr(mod, factory)(0, 1, "model.pt", task=task)
        assert mock_create.call_args.kwargs["task_name"] == expected

    @pytest.mark.parametrize("calc_path, cls_name", [
        ("ensemble_analyzer.calculators.uma", "UMAMlCalc"),
        ("ensemble_analyzer.calculators.fairchem", "FAIRChemMlCalc"),
    ])
    def test_calculator_forwards_protocol_task(self, mock_conformer, mock_protocol, calc_path, cls_name):
        import importlib
        mock_protocol.charge, mock_protocol.mult = 0, 1
        mock_protocol.functional = "model.pt"
        mock_protocol.task = "omc"
        cls = getattr(importlib.import_module(calc_path), cls_name)
        with patch("enan_calculators.get_ase_calculator") as mock_get:
            cls(mock_protocol, 4, mock_conformer)._get_ml_calculator()
        assert mock_get.call_args.kwargs["task"] == "omc"


class TestMACETask:
    @pytest.fixture
    def model_file(self, tmp_path):
        path = tmp_path / "multihead.pt"
        path.touch()
        return path

    def _create(self, method, task):
        from enan_calculators._mace import create_mace_calc
        mock_mace = MagicMock()
        with patch("enan_calculators._mace.MACECalculator", mock_mace), \
             patch("enan_calculators._mace.torch") as mock_torch:
            mock_torch.cuda.is_available.return_value = False
            create_mace_calc(0, 1, method, task=task)
        return mock_mace

    def test_head_passed_to_model(self, model_file):
        mock_mace = self._create(str(model_file), task="omol")
        assert mock_mace.call_args.kwargs["head"] == "omol"

    def test_no_head_when_task_unset(self, model_file):
        mock_mace = self._create(str(model_file), task=None)
        assert "head" not in mock_mace.call_args.kwargs

    def test_heads_are_cached_separately(self, model_file):
        from enan_calculators import _mace
        with patch.dict(_mace._PREDICTOR_CACHE, clear=True):
            self._create(str(model_file), task="omol")
            self._create(str(model_file), task="omat_pbe")
            heads = {key[-1] for key in _mace._PREDICTOR_CACHE}
        assert heads == {"omol", "omat_pbe"}

    def test_head_passed_to_foundation(self):
        from enan_calculators import _mace
        with patch.dict(_mace._PREDICTOR_CACHE, clear=True), \
             patch("enan_calculators._mace._load_foundation") as mock_load:
            self._create("mp", task="omat_pbe")
        assert mock_load.call_args.kwargs == {"head": "omat_pbe"}

    def test_calculator_forwards_protocol_task(self, mock_conformer, mock_protocol):
        from ensemble_analyzer.calculators.mace import MACEMlCalc
        mock_protocol.charge, mock_protocol.mult = 0, 1
        mock_protocol.functional = "multihead.pt"
        mock_protocol.task = "omol"
        with patch("enan_calculators.get_ase_calculator") as mock_get:
            MACEMlCalc(mock_protocol, 4, mock_conformer)._get_ml_calculator()
        assert mock_get.call_args.kwargs["task"] == "omol"


class TestAIMNetTask:
    def test_factory_raises_with_task(self):
        from enan_calculators._aimnet import create_aimnet_calc
        with pytest.raises(ValueError, match="no task/head"):
            create_aimnet_calc(0, 1, "aimnet2", task="omol")

    def test_calculator_raises_with_task(self, mock_conformer, mock_protocol):
        from ensemble_analyzer.calculators.aimnet import AIMNetCalc
        mock_protocol.charge, mock_protocol.mult = 0, 1
        mock_protocol.functional = "aimnet2"
        mock_protocol.task = "omol"
        with pytest.raises(ValueError, match="no task/head"):
            AIMNetCalc(mock_protocol, 4, mock_conformer)._get_ml_calculator()

    def test_factory_without_task(self, tmp_path):
        from enan_calculators import _aimnet
        model = tmp_path / "aimnet2.pt"
        model.touch()
        mock_ase = MagicMock()
        with patch.dict(_aimnet._PREDICTOR_CACHE, clear=True), \
             patch("enan_calculators._aimnet.AIMNet2ASE", mock_ase), \
             patch("enan_calculators._aimnet.torch"):
            _aimnet.create_aimnet_calc(0, 1, str(model))
        mock_ase.assert_called_once_with(str(model), charge=0, mult=1)
