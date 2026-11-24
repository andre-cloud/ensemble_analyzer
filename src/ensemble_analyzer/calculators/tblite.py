from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("tblite")
class TBLiteCalc(BaseMlCalc):
    label = "tblite"
    default_method = "GFN2-xTB"
    uses_task = False
