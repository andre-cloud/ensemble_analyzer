from ensemble_analyzer.calculators._ml_base import BaseMlCalc
from ensemble_analyzer.calculators.base import register_calculator


@register_calculator("skala")
class SkalaCalc(BaseMlCalc):
    label = "skala"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        return get_ase_calculator(
            "skala",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=self.protocol.functional,
            basis=self.protocol.basis,
        )
