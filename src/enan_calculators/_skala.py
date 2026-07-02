import os
from pathlib import Path
from enan_calculators._models import get_models_dir

import torch
if hasattr(torch.serialization, "add_safe_globals"):
    torch.serialization.add_safe_globals([slice])

try:
    from skala.ase import Skala
    from skala.functional import load_functional
    skala_model = True
except ImportError as e:
    print(f"REAL IMPORT ERROR: {e}")
    skala_model = None
    load_functional = None
    raise


_PREDICTOR_CACHE = {}


def create_skala_calc(charge, mult, method, basis, solvent=None):
    if skala_model is None:
        raise ImportError(
            "skala module missing. Install via: pip install skala"
        )

    model_path = get_models_dir("skala", create=False) / method
    print(f'{model_path = }')
    if not model_path.exists():
        model_path_fun = Path(str(model_path) + ".fun")
        if model_path_fun.exists():
            model_path = model_path_fun
        else:
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(
                    f"Skala model not found: {model_path}. "
                    f"Please place the downloaded weights in {get_models_dir('skala', create=False)}."
                )

    checkpoint_file = load_functional(str(model_path))

    cache_key = (method, basis, charge, mult)
    if cache_key not in _PREDICTOR_CACHE:
        _PREDICTOR_CACHE[cache_key] = Skala(
            model=checkpoint_file, xc=method.strip('.fun'), basis=basis, charge=charge, multiplicity=mult,
        )

    return _PREDICTOR_CACHE[cache_key]