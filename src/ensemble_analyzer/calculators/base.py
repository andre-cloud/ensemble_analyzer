from __future__ import annotations

from abc import ABC, abstractmethod

from typing import Callable, Tuple, Any, Optional
import importlib
import numpy as np
import os
import pkgutil


class LazyCalculatorRegistry:
    """Dict-like registry that imports calculator modules on demand.

    Convention: calculator name ``"foo"`` → module ``calculators/foo.py``.
    The module's ``@register_calculator`` decorator populates the registry at
    import time; importing is deferred until the calculator is actually needed.
    """

    def __init__(self) -> None:
        self._loaded: dict[str, type] = {}
        self._all_names: set[str] | None = None

    # -- discovery ---------------------------------------------------------------

    def _discover_all(self) -> None:
        if self._all_names is not None:
            return
        pkg = importlib.import_module("ensemble_analyzer.calculators")
        self._all_names = {
            m.name
            for m in pkgutil.iter_modules(pkg.__path__)
            if m.name not in ("base", "_ml_base", "__init__")
        }

    def _load(self, name: str) -> None:
        name = name.lower()
        if name not in self._loaded:
            try:
                importlib.import_module(f"ensemble_analyzer.calculators.{name}")
            except ImportError: 
                importlib.import_module(f"ensemble_analyzer.calculators.{name}_calc")


    # -- Mapping interface -------------------------------------------------------

    def __getitem__(self, name: str) -> type:
        name = name.lower()
        self._load(name)
        if name not in self._loaded:
            raise KeyError(
                f"Calculator {name!r} not found or not registered"
            )
        return self._loaded[name]

    def __setitem__(self, name: str, cls: type) -> None:
        self._loaded[name.lower()] = cls

    def __contains__(self, name: str) -> bool:
        name = name.lower()
        if name in self._loaded:
            return True
        try:
            self._load(name)
            return name in self._loaded
        except ImportError:
            return False

    def __iter__(self):
        self._discover_all()
        return iter(self._all_names)

    def __len__(self):
        self._discover_all()
        return len(self._all_names)

    def keys(self):
        self._discover_all()
        return iter(self._all_names)

    def get(self, name: str, default: Any = None) -> Any:
        try:
            return self[name]
        except KeyError:
            return default


def register_calculator(name: str) -> Callable:
    def decorator(cls: type) -> type:
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS: set[str] = {"tblite", "aimnet", "uma", "fairchem", "mace", "skala"}


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
            if str(i) in self.conf.energies:
                er = self.conf.energies[str(i)]
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


CALCULATOR_REGISTRY = LazyCalculatorRegistry()
