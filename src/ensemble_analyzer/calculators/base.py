from __future__ import annotations

from typing import Callable, Tuple, Any, Optional
import importlib
import os


class _Registry(dict):
    """Calculator name -> class, importing ``calculators/<name>.py`` on first use.

    Each module fills the registry through its ``@register_calculator``.
    """

    def __missing__(self, name: str) -> type:
        if name in CALCULATORS:
            importlib.import_module(f"{__package__}.{name}")
        if not dict.__contains__(self, name):
            raise KeyError(f"Calculator {name!r} not found or not registered")
        return dict.get(self, name)

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and name.lower() in CALCULATORS

    def __getitem__(self, name: str) -> type:
        return super().__getitem__(name.lower())

    def keys(self):
        return iter(CALCULATORS)


def register_calculator(name: str) -> Callable:
    def decorator(cls: type) -> type:
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls

    return decorator


ML_CALCULATORS: set[str] = {"tblite", "aimnet", "uma", "fairchem", "mace", "skala"}
CALCULATORS: tuple[str, ...] = ("orca", "gaussian", "nwchem", *sorted(ML_CALCULATORS))


class BaseCalc:
    def __init__(self, protocol, cpu: int, conf=None, temperature=298.15,
                 linear=False, cut_off=100, alpha=4, P=101.325):
        self.protocol = protocol
        self.cpu = cpu
        self.conf = conf
        self.constraints = protocol.constraints
        self.temperature = temperature
        self.linear = linear
        self.cut_off = cut_off
        self.alpha = alpha
        self.P = P

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

    def _add_nmr_keywords(self, calc: Any) -> None:
        pass

    def _nmr_settings(self):
        """The step's NMR settings, or None (protocols in tests may be MagicMocks)."""
        nmr = getattr(self.protocol, "nmr", None)
        return nmr if isinstance(nmr, dict) and nmr else None

    # --- Template methods ---

    def single_point(self) -> Tuple[Any, str]:
        calc, label = self._build_calculator()
        self._add_sp_keywords(calc)
        self._add_tddft_keywords(calc)
        self._add_nmr_keywords(calc)
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


CALCULATOR_REGISTRY = _Registry()
