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
    "omol": "mace_omol",
}


def _uses_charge_spin(calc) -> bool:
    """True if the loaded MACE model reads atoms.info["charge"/"spin"].

    MACE-OMOL-style models declare total_charge/total_spin embeddings in
    ``embedding_specs``; MACE-POLAR models (PolarMACE) take the total charge
    and spin in their electrostatics. Others (MACE-OFF, MACE-MP, Egret-1, ...)
    ignore them.
    """
    model = calc.models[0]
    if type(model).__name__ == "PolarMACE":
        return True
    specs = getattr(model, "embedding_specs", None)
    return isinstance(specs, dict) and "total_charge" in specs


def _is_parameter_free_gate_key(key: str, state: dict) -> bool:
    """True for the e3nn ``Gate`` tensors that newer MACE gate blocks drop.

    Models saved with e3nn 0.4.4 (e.g. MACE-POLAR-1) use ``e3nn.nn.Gate``,
    whose ``mul`` is an ElementwiseTensorProduct with an empty ``weight`` and
    an all-ones ``output_mask``. mace-torch >= 0.3.16 rebuilds that layer as a
    parameter-free ``GatedEquivariantBlock`` computing the same function
    (checked to 1e-15 on MACE-POLAR-1 S/M/L), so nothing learned is lost.
    """
    if key.endswith("equivariant_nonlin.mul.weight"):
        return state[key].numel() == 0
    if key.endswith("equivariant_nonlin.mul.output_mask"):
        return bool((state[key] == 1).all())
    return False


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
        if isinstance(model, torch.jit.ScriptModule):
            # compiled (TorchScript) models, e.g. Egret-1, do not unpickle
            # e3nn modules and need no rebuild
            rebuilt.append(model)
            continue
        tensor = next(iter(model.parameters()), None)
        if tensor is None:
            tensor = next(iter(model.buffers()), None)
        device = tensor.device if tensor is not None else torch.device("cpu")
        new = model.__class__(**extract_config_mace_model(model))
        state = model.state_dict()
        missing, unexpected = new.load_state_dict(state, strict=False)
        unexplained = [k for k in missing if not k.endswith("_zeroed")]
        unexpected = [k for k in unexpected if not _is_parameter_free_gate_key(k, state)]
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


def _load_model_file(model_path, device, default_dtype, head_kwargs) -> Calculator:
    """Load a local MACE model file, picking the calculator type from the model.

    Electrostatic MACE-POLAR models (mace-torch >= 0.3.16) need
    ``model_type="PolarMACE"``; they are recognised from the class of the
    loaded model, so the protocol only names the file.
    """
    calc = MACECalculator(
        model_paths=str(model_path), device=device,
        default_dtype=default_dtype, **head_kwargs,
    )
    if type(calc.models[0]).__name__ == "PolarMACE":
        calc = MACECalculator(
            model_paths=str(model_path), device=device,
            default_dtype=default_dtype, model_type="PolarMACE", **head_kwargs,
        )
    return calc


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
            loaded = _load_foundation(method, device, default_dtype, **head_kwargs)
        else:
            loaded = _load_model_file(model_path, device, default_dtype, head_kwargs)
        if _NEEDS_REBUILD:
            _rebuild_models_for_e3nn(loaded, default_dtype)
        # cache only after a successful rebuild: a failed one must not leave
        # an unusable model behind for the next call
        _PREDICTOR_CACHE[cache_key] = loaded

    calc = _PREDICTOR_CACHE[cache_key]

    if (charge != 0 or mult != 1) and not _uses_charge_spin(calc):
        import warnings
        warnings.warn(
            f"MACE model '{method}' ignores charge and spin: charge={charge}, "
            f"mult={mult} will be computed as a neutral closed-shell system. "
            "Use a charge/spin-aware model (MACE-OMOL, MACE-POLAR) for charged "
            "or open-shell species.",
            stacklevel=2,
        )

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
