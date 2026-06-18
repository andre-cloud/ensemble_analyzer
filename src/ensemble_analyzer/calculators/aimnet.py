from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("aimnet")
class AIMNetCalc(BaseMlCalc):
    label = "aimnet"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional or "aimnet2")
        return get_ase_calculator(
            "aimnet",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
        )
