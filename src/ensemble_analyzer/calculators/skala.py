from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("skala")
class SkalaCalc(BaseMlCalc):
    label = "skala"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional)
        return get_ase_calculator(
            "skala",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            basis=self.protocol.basis,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
        )
