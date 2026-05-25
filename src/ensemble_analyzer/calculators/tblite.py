from typing import Tuple, Any
from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("tblite")
class TBLiteCalc(BaseMlCalc):
    """
    Calculator wrapper for the TBLite semi-empirical method (GFN-xTB family).
    """

    label = "tblite"

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        """
        Build and return the TBLite ASE calculator.

        Args:
            **kwargs: Additional keyword arguments forwarded to the TBLite
                constructor.

        Returns:
            Any: TBLite ASE calculator instance.
        """
        try:
            from tblite.ase import TBLite as TB
        except ImportError:
            raise ImportError("tblite module missing. Install via: pip install tblite")

        method = kwargs.pop("method", self.protocol.functional or "GFN2-xTB")
        solv = None
        if self.protocol.solvent and self.protocol.solvent.solvent:
            solv = ("alpb", self.protocol.solvent.solvent)

        return TB(
            method=method,
            charge=self.protocol.charge,
            multiplicity=self.protocol.mult,
            solvation=solv,
            **kwargs,
        )
