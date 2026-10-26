from enan_calculators._ml_inference import create_ml_calc
from enan_calculators._uma import DEFAULT_TASK


def create_fairchem_calc(charge, mult, method, solvent=None, task=None):
    return create_ml_calc(
        "fairchem", charge, mult, method,
        task_name=task or DEFAULT_TASK, solvent=solvent,
    )


from enan_calculators._uma import UMAMlCalc
from ensemble_analyzer.calculators.base import register_calculator


@register_calculator("fairchem")
class FAIRChemMlCalc(UMAMlCalc):
    label = "fairchem"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional)
        if method is None:
            raise ValueError("fairchem calculator requires a model path via protocol.functional or method= kwarg")
        return get_ase_calculator(
            "fairchem",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
            task=self.protocol.task,
        )
