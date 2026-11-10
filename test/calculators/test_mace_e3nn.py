"""
MACE models saved with e3nn 0.4.4 loaded under e3nn >= 0.5 (enan_calculators._mace).

- ``_normalise_codegen_state``: e3nn 0.4.4 stores compiled submodules as raw
  TorchScript bytes, e3nn >= 0.5 expects (buffer_type, buffer) tuples.
- ``_rebuild_models_for_e3nn``: models are re-instantiated from their own
  hyper-parameters and the weights copied back.
"""

from types import SimpleNamespace

import pytest

from enan_calculators._mace import _normalise_codegen_state


class TestNormaliseCodegenState:

    def test_raw_bytes_become_torchscript_tuples(self):
        state = {"_compiled_main": b"\x00ts", "other": bytearray(b"ab")}
        _normalise_codegen_state(state)
        assert state == {"_compiled_main": ("torchscript", b"\x00ts"), "other": ("torchscript", b"ab")}

    def test_long_tuples_truncated(self):
        state = {"f": ("fx", b"buf", "extra")}
        _normalise_codegen_state(state)
        assert state == {"f": ("fx", b"buf")}

    def test_current_format_untouched(self):
        state = {"f": ("torchscript", b"buf")}
        _normalise_codegen_state(state)
        assert state == {"f": ("torchscript", b"buf")}


class TestRebuildModels:

    @pytest.fixture
    def toy_model(self, monkeypatch):
        torch = pytest.importorskip("torch")
        scripts_utils = pytest.importorskip("mace.tools.scripts_utils")
        monkeypatch.setattr(scripts_utils, "extract_config_mace_model", lambda model: {"n": 2})

        class Toy(torch.nn.Module):
            def __init__(self, n=2):
                super().__init__()
                self.w = torch.nn.Parameter(torch.ones(n, dtype=torch.float32))
                self.register_buffer("mask_zeroed", torch.zeros(n))

        old = Toy()
        with torch.no_grad():
            old.w.fill_(3.0)
        return torch, old

    def test_weights_copied_and_frozen(self, toy_model):
        from enan_calculators._mace import _rebuild_models_for_e3nn
        torch, old = toy_model
        del old._buffers["mask_zeroed"]  # buffer recomputed by the constructor, absent in old checkpoints
        calc = SimpleNamespace(models=[old])
        _rebuild_models_for_e3nn(calc, "float64")
        (new,) = calc.models
        assert new is not old
        assert torch.allclose(new.w, torch.full((2,), 3.0, dtype=torch.float64))
        assert new.w.dtype == torch.float64
        assert not new.w.requires_grad
        assert not new.training

    def test_unexplained_missing_weights_raise(self, toy_model):
        from enan_calculators._mace import _rebuild_models_for_e3nn
        torch, old = toy_model
        del old._parameters["w"]
        with pytest.raises(RuntimeError, match="could not be rebuilt"):
            _rebuild_models_for_e3nn(SimpleNamespace(models=[old]), "float64")
