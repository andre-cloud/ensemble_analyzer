from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("uma")
class UMAMlCalc(BaseMlCalc):
    label = "uma"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional or "uma-s-1.pt")
        return get_ase_calculator(
            "uma",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
        )
