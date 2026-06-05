from abc import ABC, abstractmethod

from typing import Callable, Dict, Tuple, Any, Optional
import numpy as np
import os


def register_calculator(name: str) -> Callable:
    def decorator(cls: type) -> type:
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS = {"tblite", "aimnet", "uma"}


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

    def _find_hessian_source_protocol(self) -> Optional[int]:
        if not self.protocol.ts:
            return None
        for i in range(self.protocol.number - 1, -1, -1):
            if i in self.conf.energies:
                er = self.conf.energies[i]
                if er.Freq is not None and len(er.Freq) > 0:
                    if not er.calculator or er.calculator == self.protocol.calculator:
                        return i
        return None

    # --- Hook methods (override in subclasses) ---

    def _build_calculator(self) -> Tuple[Any, str]:
        raise NotImplementedError

    def _add_sp_keywords(self, calc: Any) -> None:
        pass

    def _add_opt_keywords(self, calc: Any) -> None:
        pass

    def _add_freq_keywords(self, calc: Any) -> None:
        pass

    def _add_tddft_keywords(self, calc: Any) -> None:
        pass

    # --- Template methods ---

    def single_point(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_sp_keywords(calc)
        self._add_tddft_keywords(calc)
        return calc, label

    def optimisation(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_opt_keywords(calc)
        self._add_tddft_keywords(calc)
        return calc, label

    def frequency(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_freq_keywords(calc)
        self._add_tddft_keywords(calc)
        return calc, label

    @staticmethod
    def _build_path(*parts: str) -> str:
        return os.path.abspath(os.path.join(*parts)).replace('\\', '/')


CALCULATOR_REGISTRY : Dict[str, BaseCalc] = {}
