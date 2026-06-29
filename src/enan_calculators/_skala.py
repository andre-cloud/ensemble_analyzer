try:
    from skala.ase import Skala as _Skala
except ImportError:
    _Skala = None

_PREDICTOR_CACHE = {}


def create_skala_calc(charge, mult, method, basis, solvent=None):
    if _Skala is None:
        raise ImportError(
            "skala module missing. Install via: pip install skala"
        )

    cache_key = (method, basis, charge, mult)
    if cache_key not in _PREDICTOR_CACHE:
        _PREDICTOR_CACHE[cache_key] = _Skala(
            xc=method, basis=basis, charge=charge, multiplicity=mult,
        )

    return _PREDICTOR_CACHE[cache_key]
