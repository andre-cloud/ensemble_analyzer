from dataclasses import dataclass
from typing import List, Optional, Union, Literal
import numpy as np
from datetime import datetime

from ensemble_analyzer._spectral.graph_default import GraphDefault

from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.protocol import Protocol
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer.constants import CHIRALS


# ---------------------------------------------------------------------------
# Lazy numba — compiled on first call, pure-numpy vectorised fallback
# ---------------------------------------------------------------------------

_NUMBA_COMPILED = {}

def _gaussian_numpy(X, x0, I, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    norm = 1.0 / (sigma * np.sqrt(2 * np.pi))
    dx = (X[:, np.newaxis] - x0[np.newaxis, :]) / sigma
    return norm * np.exp(-0.5 * dx ** 2) @ I


def _lorentzian_numpy(X, x0, I, fwhm):
    n_peaks = x0.shape[0]
    n_x = X.shape[0]
    Y = np.zeros(n_x)
    if n_peaks == 0:
        return Y
    fwhm2 = fwhm * fwhm
    for i in range(n_peaks):
        xi = x0[i]
        Ii = I[i]
        if Ii == 0.0:
            continue
        dx = X - xi
        Y += Ii * fwhm2 / (fwhm2 + 4.0 * dx * dx)
    return Y


def _diversity_numpy(a, b, weight, max_val):
    diff = a - b
    return np.sqrt(np.sum(weight * diff * diff) / diff.shape[0]) / max_val


def _compile_numba():
    if _NUMBA_COMPILED:
        return
    from numba import njit, prange

    @njit(cache=True, fastmath=True)
    def _gauss(X, x0, I, fwhm):
        n_x = X.shape[0]
        n_peaks = x0.shape[0]
        Y = np.zeros(n_x)
        if n_peaks == 0:
            return Y
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        norm = 1.0 / (sigma * np.sqrt(2 * np.pi))
        inv_sigma = 1.0 / sigma
        for j in prange(n_x):
            yj = 0.0
            Xj = X[j]
            for i in range(n_peaks):
                Ii = I[i]
                if Ii == 0.0:
                    continue
                dx = (Xj - x0[i]) * inv_sigma
                yj += Ii * norm * np.exp(-0.5 * dx * dx)
            Y[j] = yj
        return Y

    @njit(cache=True, fastmath=True)
    def _lor(X, x0, I, fwhm):
        n_peaks = x0.shape[0]
        n_x = X.shape[0]
        Y = np.zeros(n_x)
        if n_peaks == 0:
            return Y
        fwhm2 = fwhm * fwhm
        for i in prange(n_peaks):
            xi = x0[i]
            Ii = I[i]
            if Ii == 0.0:
                continue
            for j in range(n_x):
                dx = X[j] - xi
                Y[j] += Ii * fwhm2 / (fwhm2 + 4.0 * dx * dx)
        return Y

    @njit(cache=True, fastmath=True)
    def _div(a, b, weight, max_val):
        diff = a - b
        s = 0.0
        n = diff.shape[0]
        for i in prange(n):
            s += weight[i] * diff[i] * diff[i]
        return np.sqrt(s / n) / max_val

    _NUMBA_COMPILED['gaussian'] = _gauss
    _NUMBA_COMPILED['lorentzian'] = _lor
    _NUMBA_COMPILED['diversity'] = _div


def gaussian_njit(X, x0, I, fwhm):
    if 'gaussian' not in _NUMBA_COMPILED:
        _compile_numba()
    return _NUMBA_COMPILED['gaussian'](X, x0, I, fwhm)


def lorentzian_njit(X, x0, I, fwhm):
    if 'lorentzian' not in _NUMBA_COMPILED:
        _compile_numba()
    return _NUMBA_COMPILED['lorentzian'](X, x0, I, fwhm)


def diversity_function_njit(a, b, weight, max_val):
    if 'diversity' not in _NUMBA_COMPILED:
        _compile_numba()
    return _NUMBA_COMPILED['diversity'](a, b, weight, max_val)


# ---------------------------------------------------------------------------
# BaseGraph
# ---------------------------------------------------------------------------


@dataclass
class BaseGraph:
    """
    Base class for spectral graph generation and convolution.

    Handles the retrieval of discrete transitions from conformers,
    convolution with line-shape functions, and auto-optimization of
    spectral parameters (shift, fwhm) against a reference.
    """

    confs: List[Conformer]
    protocol: Protocol
    graph_type: Literal['IR', 'VCD', 'UV', 'ECD']
    log: Logger

    invert: Optional[bool] = False
    fwhm_user: Optional[Union[List[float], float]] = None
    shift_user: Optional[Union[List[float], float]] = None

    read_population: Optional[Union[int, float, str]] = None
    definition: Optional[int] = 4
    interested_area: Optional[list] = None

    def __post_init__(self) -> None:
        self.defaults = GraphDefault(self.graph_type)

        self.X = np.linspace(self.defaults.start, self.defaults.end, num=10 ** self.definition)
        self.X = self.X[np.argsort(self.X)]

    def retrieve_data(self, protocol: Protocol) -> None:
        self.impulse = []
        self.energies = []
        population_from = str(self.read_population) if self.read_population else str(protocol.number)

        for conf in self.confs:
            if not self.check_conf(conf, protocol):
                continue

            p = conf.energies[population_from].Pop
            x = np.array(conf.graphs_data[protocol.number, self.graph_type].X)
            y = np.array(conf.graphs_data[protocol.number, self.graph_type].Y) * p

            if x.size < 1:
                continue
            if self.invert:
                y *= -1

            self.energies.append(x)
            self.impulse.append(y)

        if len(self.energies) > 0:
            self.energies = np.concatenate(self.energies)
            self.impulse = np.concatenate(self.impulse)
        else:
            self.energies = np.array([])
            self.impulse = np.array([])

    def normalize(self, Y: np.ndarray, idx_min: Optional[int] = None,
                  idx_max: Optional[int] = None) -> np.ndarray:
        if idx_min is not None and idx_max is not None:
            max_value = np.max(np.abs(Y[idx_min:idx_max]))
        else:
            max_value = np.max(np.abs(Y))

        return Y / max_value

    def dump_XY_data(self, X: np.ndarray, Y: np.ndarray, fname: str) -> None:
        data = np.column_stack((X, Y))
        np.savetxt(fname, data, delimiter=' ')

    def check_conf(self, conf: Conformer, protocol: Protocol) -> bool:
        if not conf.active:
            return False
        if protocol.number not in conf.graphs_data:
            return False
        if not conf.graphs_data.has_graph_type(protocol.number, self.graph_type):
            return False
        return True

    def diversity_function(self, a: np.ndarray, b: np.ndarray,
                           w: Optional[np.ndarray] = None) -> float:
        MAX = 1 if self.graph_type not in CHIRALS else 2
        w = self.ref.weight if w is None else w
        return diversity_function_njit(a=a, b=b, weight=w, max_val=MAX)

    def set_boundaries(self) -> None:
        if isinstance(self.shift_user, list):
            self.shift_bounds = self.shift_user
        elif isinstance(self.shift_user, float) or isinstance(self.shift_user, int):
            self.shift_bounds = [self.shift_user, self.shift_user]
        elif not self.shift_user:
            self.shift_bounds = self.defaults.shift_intervals

        if isinstance(self.fwhm_user, list):
            self.fwhm_bounds = self.fwhm_user
        elif isinstance(self.fwhm_user, float) or isinstance(self.fwhm_user, int):
            self.fwhm_bounds = [self.fwhm_user, self.fwhm_user]
        elif not self.fwhm_user:
            self.fwhm_bounds = self.defaults.fwhm_intervals

    def compute_spectrum(self) -> None:
        self.log.debug("Compute spectrum")
        self.set_boundaries()
        self.log.debug("Retrieving data")
        self.retrieve_data(self.protocol)

        if self.energies.size == 0 or self.energies[self.energies != 0].size == 0:
            self.log.spectra_skip(self.graph_type)
            return

        if self.ref:
            self.autoconvolution()
        else:
            self.SHIFT = self.defaults.shift
            self.FWHM = self.defaults.fwhm

            Y = self.convolute(energies=self.energies, impulses=self.impulse,
                               shift=self.SHIFT, fwhm=self.FWHM)

            self.Y = self.normalize(Y)

            self.log.spectra_result(
                graph_type=self.graph_type,
                parameters={"Shift": self.SHIFT, "FWHM": self.FWHM},
                msg=f"Using default parameters, Reference {self.graph_type} Spectra not found",
            )

        if self.Y[~np.isnan(self.Y)].size > 0:
            self.log.debug(f'Saving {self.graph_type} spectra convoluted')
            self.dump_XY_data(self.X, self.Y,
                              f'{self.graph_type}_p{self.protocol.number}_comp.xy')

    def autoconvolution(self) -> None:
        from scipy.optimize import minimize

        ref_norm = self.ref.Y

        def callback_optimizer(params):
            shift, fwhm = params
            Y_conv = self.convolute(self.energies, self.impulse, shift, fwhm)
            Y_conv = self.normalize(Y_conv, idx_min=self.ref.x_min_idx,
                                    idx_max=self.ref.x_max_idx)
            rmsd = self.diversity_function(Y_conv, ref_norm)
            return rmsd

        initial_guess = [
            sum(self.shift_bounds) * .5,
            sum(self.fwhm_bounds) * .5,
        ]

        st = datetime.now()

        result = minimize(
            fun=callback_optimizer, x0=initial_guess,
            bounds=(self.shift_bounds, self.fwhm_bounds),
            options={"maxiter": 1000},
        )
        end = datetime.now()

        if result.success:
            self.SHIFT, self.FWHM = result.x
            t = "Spectra convolution results:"
        else:
            self.SHIFT, self.FWHM = self.defaults.shift, self.defaults.fwhm
            t = "Spectra convolution did NOT converged. Using default parameters:"

        Y = self.convolute(energies=self.energies, impulses=self.impulse,
                           shift=self.SHIFT, fwhm=self.FWHM)
        self.Y = self.normalize(Y, idx_min=self.ref.x_min_idx,
                                idx_max=self.ref.x_max_idx)

        diversity = self.diversity_function(
            self.Y[self.ref.x_min_idx:self.ref.x_max_idx],
            ref_norm[self.ref.x_min_idx:self.ref.x_max_idx],
        )
        similarity = ((1 if self.graph_type not in CHIRALS else 2) - diversity) / \
                     (1 if self.graph_type not in CHIRALS else 2) * 100

        diversity_unw = self.diversity_function(
            self.Y[self.ref.x_min_idx:self.ref.x_max_idx],
            ref_norm[self.ref.x_min_idx:self.ref.x_max_idx],
            w=np.ones_like(self.Y[self.ref.x_min_idx:self.ref.x_max_idx]),
        )
        similarity_unw = ((1 if self.graph_type not in CHIRALS else 2) - diversity_unw) / \
                          (1 if self.graph_type not in CHIRALS else 2) * 100

        self.log.spectra_result(
            graph_type=self.graph_type,
            parameters={
                "Shift": self.SHIFT, "FWHM": self.FWHM,
                "Similarity": similarity,
                "Similarity Unweighted": similarity_unw,
                "Time": (end - st), "Cycle": f"{result.nfev}",
            },
            msg=t,
        )

    def gaussian(self, x0: np.ndarray, I: np.ndarray, fwhm: float) -> np.ndarray:
        return gaussian_njit(self.X, x0, I, fwhm)

    def lorentzian(self, x0: np.ndarray, I: np.ndarray, fwhm: float) -> np.ndarray:
        return lorentzian_njit(self.X, x0, I, fwhm)
