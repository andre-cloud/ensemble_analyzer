import os
from pathlib import Path
from enan_calculators._models import get_models_dir

import torch
if hasattr(torch.serialization, "add_safe_globals"):
    torch.serialization.add_safe_globals([slice])

try:
    from skala.ase import Skala
    skala_model = True
except ImportError as e:
    print(f"REAL IMPORT ERROR: {e}")
    skala_model = None
    raise


_PREDICTOR_CACHE = {}


def create_skala_calc(charge, mult, method, basis, solvent=None):
    if skala_model is None:
        raise ImportError(
            "skala module missing. Install via: pip install skala"
        )

    model_path = get_models_dir("skala", create=False) / method
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

    cache_key = (method, basis, charge, mult)

    if cache_key not in _PREDICTOR_CACHE:
        os.environ["SKALA_LOCAL_MODEL_PATH"] = str(model_path)

        c = Skala(
            ks_config={"functional_path": model_path}, 
            xc=str(method).removesuffix('.fun'), 
            basis=basis, 
            charge=charge, 
            multiplicity=mult,
            with_density_fit=True,
        )

        _PREDICTOR_CACHE[cache_key] = c

    return _PREDICTOR_CACHE[cache_key]