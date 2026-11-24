from enan_calculators._ml_inference import create_ml_calc
from enan_calculators._uma import DEFAULT_TASK


def create_fairchem_calc(charge, mult, method, solvent=None, task=None):
    return create_ml_calc("fairchem", charge, mult, method, task_name=task or DEFAULT_TASK)
