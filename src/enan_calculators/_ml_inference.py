from pathlib import Path
from enan_calculators._models import get_models_dir

try:
    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
    import torch
except ImportError:
    FAIRChemCalculator = None
    load_predict_unit = None
    InferenceSettings = None
    torch = None

_PREDICTOR_CACHE = {}
_FAIRCHEM_DIRS = ("fairchem", "uma")


def _get_cpu_threads():
    import os
    return (
        int(os.environ.get("SLURM_CPUS_PER_TASK", 0))
        or int(os.environ.get("OMP_NUM_THREADS", 0))
        or os.cpu_count()
    )


def create_ml_calc(calc_name, charge, mult, method, task_name, merge_mole=False, solvent=None):
    if FAIRChemCalculator is None:
        raise ImportError("fairchem-core missing. Install via: pip install \"ensemble-analyzer[fairchem]\"")

    if hasattr(torch.serialization, "add_safe_globals"):
        torch.serialization.add_safe_globals([slice])

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        torch.set_float32_matmul_precision("high")
    else:
        torch.set_num_threads(_get_cpu_threads())

    # "uma" and "fairchem" are the same fairchem-core (v2) backend: look in
    # both model folders so weights work whichever engine name is used.
    search_dirs = [get_models_dir(name, create=False) for name in (calc_name, *_FAIRCHEM_DIRS)]
    candidates = [d / method for d in dict.fromkeys(search_dirs)] + [Path(method)]
    model_path = next((c for c in candidates if c.exists()), None)
    if model_path is None:
        raise FileNotFoundError(
            f"Model file not found: {method}. "
            f"Please place the downloaded weights in {get_models_dir('fairchem', create=False)}."
        )

    # Merged MoLE weights are valid for one composition only, so a merged
    # predictor is never shared between different systems.
    cache_key = (str(model_path.resolve()), device, merge_mole)

    if cache_key not in _PREDICTOR_CACHE:
        inference_settings = InferenceSettings(
            tf32=True,
            activation_checkpointing=False,
            merge_mole=merge_mole,
            compile=True,
            max_atoms=256,
        )

        print(f"Loading and compiling {calc_name} predictor ({method}) on {device}...")
        _PREDICTOR_CACHE[cache_key] = load_predict_unit(
            path=model_path,
            device=device,
            inference_settings=inference_settings,
        )

    predictor = _PREDICTOR_CACHE[cache_key]
    if merge_mole:
        # merging is bound to the first system seen: never reuse it
        _PREDICTOR_CACHE.pop(cache_key)

    from ase.calculators.calculator import all_changes

    class _MLWrappedCalc(FAIRChemCalculator):
        def __init__(self, predictor_unit, _charge, _mult, _task_name):
            super().__init__(predictor_unit, task_name=_task_name)
            self._ml_charge = _charge
            self._ml_mult = _mult

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            calc_atoms = atoms if atoms is not None else self.atoms
            if calc_atoms is not None:
                calc_atoms.info["charge"] = self._ml_charge
                calc_atoms.info["spin"] = self._ml_mult
            super().calculate(atoms, properties, system_changes)

    return _MLWrappedCalc(predictor, charge, mult, task_name)
