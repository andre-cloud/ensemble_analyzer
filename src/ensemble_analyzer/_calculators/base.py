from abc import ABC, abstractmethod

from typing import Dict, Tuple, Any


def register_calculator(name):
    """Decorator to register each calculator class."""

    def decorator(cls):
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


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

CALCULATOR_REGISTRY : Dict[str, BaseCalc]= {}