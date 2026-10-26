from enan_calculators._ml_inference import create_ml_calc

DEFAULT_TASK = "omol"


def create_uma_calc(charge, mult, method, solvent=None, task=None):
    return create_ml_calc(
        "uma", charge, mult, method,
        task_name=task or DEFAULT_TASK, merge_mole=True, solvent=solvent,
    )


from ensemble_analyzer.calculators._ml_base import BaseMlCalc
from ensemble_analyzer.calculators.base import register_calculator


@register_calculator("uma")
class UMAMlCalc(BaseMlCalc):
    label = "uma"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional)
        if method is None:
            raise ValueError("uma calculator requires a model path via protocol.functional or method= kwarg")
        return get_ase_calculator(
            "uma",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
            task=self.protocol.task,
        )
