from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("uma")
class UMAMlCalc(BaseMlCalc):
    label = "uma"
    method_required = "a model path via protocol.functional or method= kwarg"
