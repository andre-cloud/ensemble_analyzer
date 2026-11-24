from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("aimnet")
class AIMNetCalc(BaseMlCalc):
    label = "aimnet"
    default_method = "aimnet2"
