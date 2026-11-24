from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("skala")
class SkalaCalc(BaseMlCalc):
    label = "skala"
    uses_task = False
    uses_basis = True
