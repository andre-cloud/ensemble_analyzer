import os
from pathlib import Path
from enan_calculators._models import get_models_dir

_Skala = None
_PREDICTOR_CACHE = {}


def _load_skala():
    global _Skala
    if _Skala is not None:
        return
    try:
        import torch
        if hasattr(torch.serialization, "add_safe_globals"):
            torch.serialization.add_safe_globals([slice])
    except Exception:
        pass
    try:
        from skala.ase import Skala
        from skala.functional import load_functional
        _Skala = Skala
    except Exception as e:
        raise ImportError(
            f"skala module error: {e}\n"
            "Install via: pip install skala"
        )


def create_skala_calc(charge, mult, method, basis, solvent=None):
    _load_skala()

    model_path = get_models_dir("skala", create=False) / method
    if not model_path.exists():
        model_path_fun = Path(str(model_path) + ".fun")
        if model_path_fun.exists():
            model_path = model_path_fun
        else:
            model_path = Path(method)

    if model_path.exists():
        skala_file = str(model_path)
        checkpoint_file = load_functional(model_path)
    else: 
        raise FileExistsError(f'{model_path} does not exists.')

    cache_key = (method, basis, charge, mult)
    if cache_key not in _PREDICTOR_CACHE:
        _PREDICTOR_CACHE[cache_key] = _Skala(
            xc=checkpoint_file, basis=basis, charge=charge, multiplicity=mult,
        )

    return _PREDICTOR_CACHE[cache_key]
