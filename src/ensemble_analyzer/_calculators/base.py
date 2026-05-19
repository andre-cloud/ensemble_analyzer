from abc import ABC, abstractmethod

from typing import Callable, Dict, Tuple, Any
import numpy as np


def register_calculator(name: str) -> Callable:
    """Decorator to register each calculator class.

    Args:
        name: Calculator identifier (e.g. 'orca', 'gaussian').

    Returns:
        Decorator that registers the class in CALCULATOR_REGISTRY.
    """

    def decorator(cls: type) -> type:
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS = {"tblite", "aimnet"}


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
    def common_str(self) -> dict:
        """Generate common input keywords for the calculator.

        Returns:
            dict: Dictionary of input keywords for the calculator.
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
    Base class for ML calculator wrappers (TBLite, AIMNet).
    Overrides common_str, optimisation, and frequency for ML behaviour.

    ML calculators use OMP_NUM_THREADS to control internal threading
    (set in calculation_executor.py before each run).
    """

    def common_str(self) -> dict:
        """ML calculators have no special input keywords.

        Returns:
            dict: Empty dictionary.
        """
        return {}

    def _get_ml_calculator(self, **kwargs: Any) -> Any:
        """Override in subclass to return the ML ASE Calculator.

        Returns:
            ASE Calculator instance.
        """
        raise NotImplementedError

    def single_point(self) -> Tuple[Any, str]:
        """Run a single-point energy calculation with the ML calculator.

        Returns:
            Tuple[Any, str]: ASE Calculator instance and label.
        """
        calc = self._get_ml_calculator()
        return calc, self.label

    def optimisation(self) -> Tuple[Any, str]:
        """Run a geometry optimisation with the ML calculator using BFGS.

        Returns:
            Tuple[Any, str]: ASE Calculator instance and label.
        """
        from ase.optimize import BFGS
        calc = self._get_ml_calculator()
        atoms = self.conf.get_ase_atoms(calc)
        with BFGS(atoms) as opt:
            opt.run(fmax=0.05)
        self.conf.last_geometry = atoms.get_positions().copy()
        return calc, self.label

    def frequency(self) -> Tuple[Any, str]:
        """Frequency calculations are not supported for ML calculators.

        Raises:
            NotImplementedError: Always raised.
        """
        raise NotImplementedError(
            f"Frequency not implemented for ML calculator '{self.label}'. "
            "Use a QM calculator (ORCA/Gaussian/NWChem) for Freq steps."
        )


CALCULATOR_REGISTRY : Dict[str, BaseCalc] = {}