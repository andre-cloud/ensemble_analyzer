from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("mace")
class MACEMlCalc(BaseMlCalc):
    label = "mace"
    method_required = (
        "a model via protocol.functional or method= kwarg "
        "(a foundation name: off, omol, mp, anicc, mdp; or a .model file)"
    )
