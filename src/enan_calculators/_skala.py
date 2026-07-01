import os
from pathlib import Path
from enan_calculators._models import get_models_dir

# Global imports for skala
try:
    from skala.functional import load_functional
    import skala.ase as skala_module
except ImportError:
    load_functional = None
    skala_module = None

_Skala = None
_PREDICTOR_CACHE = {}

def _load_skala():
    global _Skala
    if _Skala is not None:
        return
    
    if skala_module is None:
        raise ImportError("skala module error: Install via: pip install skala")
        
    try:
        import torch
        if hasattr(torch.serialization, "add_safe_globals"):
            torch.serialization.add_safe_globals([slice])
    except Exception:
        pass
        
    _Skala = skala_module.Skala

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
        # load_functional is now globally available
        checkpoint_file = load_functional(model_path)
    else: 
        raise FileExistsError(f'{model_path} does not exist.')

    cache_key = (method, basis, charge, mult)
    if cache_key not in _PREDICTOR_CACHE:
        _PREDICTOR_CACHE[cache_key] = _Skala(
            xc=checkpoint_file, basis=basis, charge=charge, multiplicity=mult,
        )

    return _PREDICTOR_CACHE[cache_key]