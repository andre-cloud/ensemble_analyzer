import os
from pathlib import Path
from ._models import get_models_dir
from ase.calculators.calculator import Calculator, all_changes

try:
    import torch
    from mace.calculators import MACECalculator
except ImportError:
    torch = None
    MACECalculator = None

_FOUNDATION = {
    "mp": "mace_mp",
    "off": "mace_off",
    "anicc": "mace_anicc",
    "mdp": "mace_mdp",
}


def _load_foundation(name: str, device: str, default_dtype: str) -> Calculator:
    func_name = _FOUNDATION[name]
    try:
        from mace.calculators import foundations_models
        func = getattr(foundations_models, func_name)
    except (ImportError, AttributeError):
        import importlib
        mod = importlib.import_module("mace.calculators")
        func = getattr(mod, func_name)
    return func(device=device, default_dtype=default_dtype)


def create_mace_calc(charge, mult, method, solvent=None):
    if MACECalculator is None:
        raise ImportError(
            "mace module missing. Install via: pip install mace-torch"
        )

    os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

    device = "cuda" if torch.cuda.is_available() else "cpu"
    default_dtype = "float64"

    if method in _FOUNDATION:
        calc = _load_foundation(method, device, default_dtype)
    else:
        model_path = get_models_dir("mace", create=False) / method
        if not model_path.exists():
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(
                    f"MACE model not found: {model_path}. "
                    f"Please place the downloaded weights in {get_models_dir('mace', create=False)}."
                )

        calc = MACECalculator(
                    model_paths=str(model_path),
                    device=device,
                    default_dtype=default_dtype,
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
