from typing import Any
from pathlib import Path
from ._ml_base import BaseMlCalc
from .base import register_calculator
from ensemble_analyzer.constants import get_models_dir
from ase.calculators.calculator import all_changes

try:
    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit
    import torch
except ImportError:
    FAIRChemCalculator = None
    load_predict_unit = None
    torch = None


@register_calculator("uma")
class UMAMlCalc(BaseMlCalc):
    label = "uma"

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        if FAIRChemCalculator is None:
            raise ImportError(
                "fairchem-core module missing. Install via: pip install fairchem-core"
            )

        method = kwargs.pop("method", self.protocol.functional or "uma-s-1.pt")
        device = "cuda" if torch.cuda.is_available() else "cpu"

        model_path = get_models_dir("uma", create=False) / method
        if not model_path.exists():
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(
                    f"Model file not found: {model_path}. "
                    f"Please place the downloaded weights in {get_models_dir('uma', create=False)}."
                )

        predictor = load_predict_unit(
            path=model_path,
            device=device,
        )

        class _UMAWrappedCalc(FAIRChemCalculator):
            def __init__(self, predictor_unit, charge, mult):
                super().__init__(predictor_unit, task_name="omol")
                self._uma_charge = charge
                self._uma_mult = mult

            def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
                calc_atoms = atoms if atoms is not None else self.atoms
                if calc_atoms is not None:
                    calc_atoms.info["charge"] = self._uma_charge
                    calc_atoms.info["spin"] = self._uma_mult
                super().calculate(atoms, properties, system_changes)

        return _UMAWrappedCalc(predictor, self.protocol.charge, self.protocol.mult)
