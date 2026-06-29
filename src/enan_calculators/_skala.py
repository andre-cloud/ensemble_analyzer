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
        from skala.ase import Skala
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
        os.environ["SKALA_LOCAL_MODEL_PATH"] = str(model_path)
    elif "SKALA_LOCAL_MODEL_PATH" in os.environ:
        del os.environ["SKALA_LOCAL_MODEL_PATH"]

    cache_key = (method, basis, charge, mult)
    if cache_key not in _PREDICTOR_CACHE:
        _PREDICTOR_CACHE[cache_key] = _Skala(
            xc=method, basis=basis, charge=charge, multiplicity=mult,
        )

    return _PREDICTOR_CACHE[cache_key]
