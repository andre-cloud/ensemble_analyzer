from ensemble_analyzer.calculators._ml_base import BaseMlCalc
from ensemble_analyzer.calculators.base import register_calculator


@register_calculator("skala")
class SkalaCalc(BaseMlCalc):
    label = "skala"

    def _get_ml_calculator(self, **kwargs):
        from enan_calculators import get_ase_calculator
        method = kwargs.pop("method", self.protocol.functional)
        return get_ase_calculator(
            "skala",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=method,
            basis=self.protocol.basis,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
        )




if __name__=='__main__':
    from ase.build import molecule
    from types import SimpleNamespace

    atoms = molecule("H2O")
    calc = SkalaCalc(
        protocol=SimpleNamespace(
            functional="skala-1.1", charge=0, mult=1, basis="def2-svp",
            solvent=None, constrains=[]
        ),
        cpu=1,
    )
    ml_calc, _ = calc.single_point()
    atoms.calc = ml_calc
    print(atoms.get_potential_energy())