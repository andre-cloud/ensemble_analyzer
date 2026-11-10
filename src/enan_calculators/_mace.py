import os
from pathlib import Path
from enan_calculators._models import get_models_dir
from ase.calculators.calculator import Calculator, all_changes

try:
    import torch
    from mace.calculators import MACECalculator
except ImportError:
    torch = None
    MACECalculator = None

_PREDICTOR_CACHE = {}


# mace-torch pins e3nn==0.4.4, FAIRChem/UMA needs e3nn>=0.5. MACE models are
# pickled with e3nn 0.4.4, whose CodeGenMixin stores each compiled submodule in
# state["__codegen__"] as raw TorchScript bytes; e3nn>=0.5 expects
# (buffer_type, buffer) tuples. Normalise the state so MACE models load with
# any e3nn version (verified: MACE-OFF23 with e3nn 0.6.0 + fairchem-core 2.25).
def _normalise_codegen_state(codegen_state: dict) -> None:
    for name, value in codegen_state.items():
        if isinstance(value, (bytes, bytearray)):
            codegen_state[name] = ("torchscript", bytes(value))
        elif isinstance(value, tuple) and len(value) > 2:
            codegen_state[name] = value[:2]


_NEEDS_REBUILD = False
try:
    import e3nn
    from e3nn.util.codegen import _mixin
    _e3nn_major_minor = tuple(int(x) for x in e3nn.__version__.split(".")[:2])
    _NEEDS_REBUILD = _e3nn_major_minor >= (0, 5)
    if _NEEDS_REBUILD and hasattr(_mixin.CodeGenMixin, "__setstate__"):
        _orig_setstate = _mixin.CodeGenMixin.__setstate__

        def _patched_setstate(self, state):
            for key in ("__codegen__", "codegen_state"):
                if isinstance(state.get(key), dict):
                    state = dict(state)
                    state[key] = dict(state[key])
                    _normalise_codegen_state(state[key])
            _orig_setstate(self, state)

        _mixin.CodeGenMixin.__setstate__ = _patched_setstate
except Exception:
    pass



_FOUNDATION = {
    "mp": "mace_mp",
    "off": "mace_off",
    "anicc": "mace_anicc",
    "mdp": "mace_mdp",
}


def _rebuild_models_for_e3nn(calc: Calculator, default_dtype: str) -> None:
    """Re-instantiate the loaded MACE models with the installed e3nn.

    Unpickling a model saved with e3nn 0.4.4 under e3nn >= 0.5 restores the
    weights but not the attributes newer e3nn modules create in ``__init__``
    (e.g. ``SphericalHarmonics.sph_func``, ``Activation.paths``). The model is
    rebuilt from its own hyper-parameters (``extract_config_mace_model``, as in
    MACE's e3nn->cueq converter) and the weights are copied back. Buffers that
    the constructor recomputes (``*_zeroed`` masks) may be absent from old
    checkpoints. Validated on MACE-OFF23: |ΔE| = 2e-5 eV vs e3nn 0.4.4.
    """
    from mace.tools.scripts_utils import extract_config_mace_model

    dtype = torch.float64 if default_dtype == "float64" else torch.float32
    rebuilt = []
    for model in calc.models:
        tensor = next(iter(model.parameters()), None)
        if tensor is None:
            tensor = next(iter(model.buffers()), None)
        device = tensor.device if tensor is not None else torch.device("cpu")
        new = model.__class__(**extract_config_mace_model(model))
        missing, unexpected = new.load_state_dict(model.state_dict(), strict=False)
        unexplained = [k for k in missing if not k.endswith("_zeroed")]
        if unexplained or unexpected:
            raise RuntimeError(
                "MACE model could not be rebuilt for e3nn "
                f"{e3nn.__version__}: missing {unexplained}, unexpected {list(unexpected)}"
            )
        new = new.to(device=device, dtype=dtype).eval()
        for param in new.parameters():
            param.requires_grad = False
        rebuilt.append(new)
    calc.models = rebuilt


def _load_foundation(name: str, device: str, default_dtype: str, **kwargs) -> Calculator:
    func_name = _FOUNDATION[name]
    try:
        from mace.calculators import foundations_models
        func = getattr(foundations_models, func_name)
    except (ImportError, AttributeError):
        import importlib
        mod = importlib.import_module("mace.calculators")
        func = getattr(mod, func_name)
    return func(device=device, default_dtype=default_dtype, **kwargs)


def create_mace_calc(charge, mult, method, solvent=None, task=None):
    if MACECalculator is None:
        raise ImportError(
            "mace module missing. Install via: pip install mace-torch"
        )

    os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

    device = "cuda" if torch.cuda.is_available() else "cpu"
    default_dtype = "float64"
    # task selects the head of a multi-head model; None keeps the model default
    head_kwargs = {"head": task} if task else {}

    if method in _FOUNDATION:
        cache_key = (method, device, default_dtype, task)
    else:
        model_path = get_models_dir("mace", create=False) / method
        if not model_path.exists():
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(
                    f"MACE model not found: {model_path}. "
                    f"Please place the downloaded weights in {get_models_dir('mace', create=False)}."
                )
        cache_key = (str(model_path), device, default_dtype, task)

    if cache_key not in _PREDICTOR_CACHE:
        if method in _FOUNDATION:
            _PREDICTOR_CACHE[cache_key] = _load_foundation(method, device, default_dtype, **head_kwargs)
        else:
            _PREDICTOR_CACHE[cache_key] = MACECalculator(
                model_paths=str(model_path),
                device=device,
                default_dtype=default_dtype,
                **head_kwargs,
            )
        if _NEEDS_REBUILD:
            _rebuild_models_for_e3nn(_PREDICTOR_CACHE[cache_key], default_dtype)

    calc = _PREDICTOR_CACHE[cache_key]

    return _MACEWrappedCalc(calc, charge, mult)


class _MACEWrappedCalc(Calculator):
    def __init__(self, calc: Calculator, _charge: int, _mult: int):
        super().__init__()
        self._inner_calc = calc
        self._mace_charge = _charge
        self._mace_mult = _mult
        self.implemented_properties = calc.implemented_properties

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        calc_atoms = atoms if atoms is not None else self.atoms
        if calc_atoms is not None:
            calc_atoms.info["charge"] = self._mace_charge
            calc_atoms.info["spin"] = self._mace_mult
        self._inner_calc.calculate(calc_atoms, properties, system_changes)
        self.results = self._inner_calc.results
