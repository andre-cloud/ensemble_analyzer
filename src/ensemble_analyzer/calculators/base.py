from abc import ABC, abstractmethod

from typing import Callable, Dict, Tuple, Any
import numpy as np
import os


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

    # --- Hook methods (override in subclasses) ---

    def _build_calculator(self) -> Tuple[Any, str]:
        raise NotImplementedError

    def _add_sp_keywords(self, calc: Any) -> None:
        pass

    def _add_opt_keywords(self, calc: Any) -> None:
        pass

    def _add_freq_keywords(self, calc: Any) -> None:
        pass

    # --- Template methods ---

    def single_point(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_sp_keywords(calc)
        return calc, label

    def optimisation(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_opt_keywords(calc)
        return calc, label

    def frequency(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_freq_keywords(calc)
        return calc, label

    @staticmethod
    def _build_path(*parts: str) -> str:
        return os.path.abspath(os.path.join(*parts)).replace('\\', '/')


CALCULATOR_REGISTRY : Dict[str, BaseCalc] = {}
