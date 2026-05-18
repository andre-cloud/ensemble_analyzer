import os
from typing import Tuple, Any
from .base import BaseMlCalc, register_calculator
from ensemble_analyzer.constants import get_models_dir


@register_calculator("aimnet")
class AIMNetCalc(BaseMlCalc):
    """
    Calculator wrapper for the AIMNet2 neural network potential.
    """

    label = "aimnet"

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        """
        Build and return the AIMNet2 ASE calculator.

        Args:
            **kwargs: Additional keyword arguments forwarded to the AIMNet2ASE
                constructor.

        Returns:
            Any: AIMNet2 ASE calculator instance.
        """
        try:
            import torch
            from aimnet.calculators import AIMNet2ASE
        except ImportError:
            raise ImportError(
                "aimnet module missing. Install via: pip install aimnet[ase]@git+https://github.com/isayevlab/aimnetcentral.git"
            )

        method = kwargs.pop("method", self.protocol.functional or "aimnet2")
        model_path = get_models_dir("aimnet", create=False) / method

        if not model_path.exists():
            model_path = Path(method)
            if not model_path.exists():
                raise FileNotFoundError(f"AIMNet model not found: {model_path}")

        torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", 1)))

        return AIMNet2ASE(
            str(model_path),
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            **kwargs,
        )
