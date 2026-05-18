from abc import ABC, abstractmethod

from typing import Dict, Tuple, Any
import numpy as np


def register_calculator(name):
    """Decorator to register each calculator class."""

    def decorator(cls):
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS = {"tblite", "aimnet", "uma"}


class BaseCalc(ABC):
    """
    Abstract Base Class for QM Calculator wrappers.
    Wraps ASE calculators to inject protocol-specific logic.
    """

    def __init__(self, protocol, cpu: int, conf=None):
        """
        Initialize the calculator wrapper.

        Args:
            protocol (Protocol): Protocol configuration object.
            cpu (int): Number of CPUs to use.
            conf (Conformer, optional): Conformer object to calculate. Defaults to None.
        """

        self.protocol = protocol
        self.cpu = cpu
        self.conf = conf
        self.constrains = protocol.constrains

    @abstractmethod
    def common_str(self):
        """
        Generate common input strings (keywords, blocks) for the calculator.

        Returns:
            Union[str, Tuple[str, str]]: Input string(s) for the calculator.
        """
        pass

    @abstractmethod
    def single_point(self) -> Tuple[Any, str]:
        """
        Configure a Single Point Energy calculation.

        Returns:
            Tuple[Calculator, str]: ASE Calculator instance and label.
        """
        pass

    @abstractmethod
    def optimisation(self) -> Tuple[Any, str]:
        """
        Configure a Geometry Optimization calculation.

        Returns:
            Tuple[Calculator, str]: ASE Calculator instance and label.
        """
        pass

    @abstractmethod
    def frequency(self) -> Tuple[Any, str]:
        """
        Configure a Frequency calculation.

        Returns:
            Tuple[Calculator, str]: ASE Calculator instance and label.
        """
        pass


class BaseMlCalc(BaseCalc):
    """
    Base class for ML calculator wrappers (TBLite, AIMNet, UMA).
    Overrides common_str, optimisation, and frequency for ML behaviour.
    """

    def common_str(self) -> str:
        return ""

    def _get_ml_calculator(self, **kwargs):
        """Override in subclass to return the ML ASE Calculator."""
        raise NotImplementedError

    def single_point(self) -> Tuple[Any, str]:
        calc = self._get_ml_calculator()
        return calc, self.label

    def optimisation(self) -> Tuple[Any, str]:
        from ase.optimize import BFGS
        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)
        with BFGS(atoms) as opt:
            opt.run(fmax=0.05)
        self.conf.last_geometry = atoms.get_positions().copy()
        return calc, self.label

    def frequency(self) -> Tuple[Any, str]:
        raise NotImplementedError(
            f"Frequency not implemented for ML calculator '{self.label}'. "
            "Use a QM calculator (ORCA/Gaussian) for freq steps."
        )


CALCULATOR_REGISTRY : Dict[str, BaseCalc] = {}