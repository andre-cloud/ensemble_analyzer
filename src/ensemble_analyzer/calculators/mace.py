from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("mace")
class MACEMlCalc(BaseMlCalc):
    label = "mace"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional)
        if method is None:
            raise ValueError(
                "mace calculator requires a model via protocol.functional or method= kwarg "
                "(a foundation name: off, omol, mp, anicc, mdp; or a .model file)"
            )
        return get_ase_calculator(
            "mace",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
            task=self.protocol.task,
        )
