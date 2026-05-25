from abc import ABC, abstractmethod

from typing import Callable, Dict, Tuple, Any
import numpy as np


def register_calculator(name: str) -> Callable:
    def decorator(cls: type) -> type:
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS = {"tblite", "aimnet"}


class BaseCalc(ABC):
    def __init__(self, protocol, cpu: int, conf=None, temperature=298.15,
                 linear=False, cut_off=100, alpha=4, P=101.325):
        self.protocol = protocol
        self.cpu = cpu
        self.conf = conf
        self.constrains = protocol.constrains
        self.temperature = temperature
        self.linear = linear
        self.cut_off = cut_off
        self.alpha = alpha
        self.P = P

    @abstractmethod
    def common_str(self) -> dict:
        pass

    @abstractmethod
    def single_point(self) -> Tuple[Any, str]:
        pass

    @abstractmethod
    def optimisation(self) -> Tuple[Any, str]:
        pass

    @abstractmethod
    def frequency(self) -> Tuple[Any, str]:
        pass


CALCULATOR_REGISTRY : Dict[str, BaseCalc] = {}
