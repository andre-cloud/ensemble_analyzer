from .base import register_calculator
from .uma import UMAMlCalc


@register_calculator("fairchem")
class FAIRChemMlCalc(UMAMlCalc):
    """Same fairchem-core backend as ``uma``, registered under its own name."""
    label = "fairchem"
