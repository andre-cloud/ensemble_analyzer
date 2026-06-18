from pathlib import Path
from ._models import get_models_dir
from ase.calculators.calculator import all_changes

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


def create_uma_calc(charge, mult, method, solvent=None):
    if FAIRChemCalculator is None:
        raise ImportError(
            "fairchem-core module missing. Install via: pip install fairchem-core"
        )

    if hasattr(torch.serialization, 'add_safe_globals'):
        torch.serialization.add_safe_globals([slice])
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        torch.set_float32_matmul_precision('high')
    else:
        import os
        torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", 1)))

    model_path = get_models_dir("uma", create=False) / method
    if not model_path.exists():
        model_path = Path(method)
        if not model_path.exists():
            raise FileNotFoundError(
                f"Model file not found: {model_path}. "
                f"Please place the downloaded weights in {get_models_dir('uma', create=False)}."
            )

    inference_settings = InferenceSettings(
        tf32=True,
        activation_checkpointing=False,
        merge_mole=True,
        compile=True,
        max_atoms=256,
    )
    predictor = load_predict_unit(
        path=model_path,
        device=device,
        inference_settings=inference_settings,
    )

    class _UMAWrappedCalc(FAIRChemCalculator):
        def __init__(self, predictor_unit, _charge, _mult):
            super().__init__(predictor_unit, task_name="omol")
            self._uma_charge = _charge
            self._uma_mult = _mult

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            calc_atoms = atoms if atoms is not None else self.atoms
            if calc_atoms is not None:
                calc_atoms.info["charge"] = self._uma_charge
                calc_atoms.info["spin"] = self._uma_mult
            super().calculate(atoms, properties, system_changes)

    return _UMAWrappedCalc(predictor, charge, mult)
