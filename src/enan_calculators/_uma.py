from enan_calculators._ml_inference import create_ml_calc

DEFAULT_TASK = "omol"


def create_uma_calc(charge, mult, method, solvent=None, task=None):
    return create_ml_calc("uma", charge, mult, method, task_name=task or DEFAULT_TASK)
