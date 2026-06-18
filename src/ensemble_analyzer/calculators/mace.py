from ._ml_base import BaseMlCalc
from .base import register_calculator


@register_calculator("mace")
class MACEMlCalc(BaseMlCalc):
    label = "mace"

    _FOUNDATION = {
        "mp": "mace_mp",
        "off": "mace_off",
        "anicc": "mace_anicc",
        "mdp": "mace_mdp",
    }

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional or "MACE_model.pt")
        return get_ase_calculator(
            "mace",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
        )
