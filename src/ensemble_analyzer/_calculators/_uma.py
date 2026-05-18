import os
from typing import Tuple, Any
from pathlib import Path
from .base import BaseMlCalc, register_calculator
from ensemble_analyzer.constants import get_models_dir
from pathlib import Path

class UMAWrappedCalc:
    """Wraps FAIRChemCalculator to inject charge/spin into atoms.info before each calculation."""

    def __init__(self, predictor_unit: Any, charge: int, mult: int) -> None:
        """
        Initialize the wrapped UMA calculator.

        Args:
            predictor_unit (Any): Loaded MLIP predictor unit.
            charge (int): Molecular charge.
            mult (int): Spin multiplicity.
        """
        from fairchem.core import FAIRChemCalculator
        self._calc = FAIRChemCalculator(predictor_unit, task_name="omol")
        self._uma_charge = charge
        self._uma_mult = mult

    def calculate(
        self,
        atoms: Any = None,
        properties: list[str] = ["energy"],
        system_changes: str = "all",
    ) -> None:
        """
        Perform calculation, injecting charge and spin into atoms info.

        Args:
            atoms (Any, optional): ASE Atoms object. Defaults to None.
            properties (list[str], optional): Properties to compute.
                Defaults to ["energy"].
            system_changes (str, optional): System changes flag.
                Defaults to "all".
        """
        calc_atoms = atoms if atoms is not None else getattr(self, 'atoms', None)
        if calc_atoms is not None:
            calc_atoms.info["charge"] = self._uma_charge
            calc_atoms.info["spin"] = self._uma_mult
        self._calc.calculate(atoms, properties, system_changes)

    def __getattr__(self, name: str) -> Any:
        """Fallback to wrapped FAIRChemCalculator attribute."""
        return getattr(self._calc, name)


@register_calculator("uma")
class UMACalc(BaseMlCalc):
    """
    Calculator wrapper for the UMA (Universal Machine-learning Atomic) potential.
    """

    label = "uma"

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        """
        Build and return the UMA ASE calculator.

        Args:
            **kwargs: Additional keyword arguments forwarded to the predictor
                loader.

        Returns:
            Any: UMA wrapped ASE calculator instance.
        """
        try:
            from fairchem.core.units.mlip_unit import load_predict_unit
            import torch
        except ImportError:
            raise ImportError("fairchem-core missing. Install via: pip install fairchem-core")

        torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", 1)))

        method = kwargs.pop("method", self.protocol.functional or "uma-s-1.pt")
        device = "cuda" if torch.cuda.is_available() else "cpu"

        model_path = get_models_dir("uma", create=False) / method

        if not model_path.exists():
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(f"UMA model not found: {model_path}")

        predictor = load_predict_unit(path=model_path, device=device, inference_settings="turbo")

        return UMAWrappedCalc(predictor, self.protocol.charge, self.protocol.mult)
