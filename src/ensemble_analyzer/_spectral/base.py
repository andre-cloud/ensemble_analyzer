from dataclasses import dataclass
from typing import List, Optional, Union, Literal
import numpy as np
import logging
from scipy.optimize import minimize
from numba import njit, prange

from datetime import datetime

from ensemble_analyzer._spectral.graph_default import GraphDefault

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer.constants import *

@dataclass
class BaseGraph: 

    confs: List[Conformer]  
    protocol : Protocol
    graph_type: Literal['IR', 'VCD', 'UV', 'ECD']
    log: Logger  

    invert : Optional[bool] = False
    fwhm_user: Optional[Union[List[float], float]] = None
    shift_user: Optional[Union[List[float], float]] = None

    read_population: Optional[Union[int, float, str]] = None
    definition: Optional[int] = 4
    interested_area: Optional[list] = None

    def __post_init__(self):
        self.defaults = GraphDefault(self.graph_type)

        self.X = np.linspace(self.defaults.start, self.defaults.end, num=10**self.definition)
        self.X = self.X[np.argsort(self.X)]

    def retrieve_data(self, protocol:Protocol) -> None: 
        self.impulse = []
        self.energies = []
        population_from = str(self.read_population) if self.read_population else str(protocol.number)

        for conf in self.confs:

            if not self.check_conf(conf, protocol):
                continue           

            p = conf.energies.__getitem__(protocol_number=population_from).Pop
            x = np.array(conf.graphs_data.__getitem__(protocol_number=protocol.number, graph_type=self.graph_type).X)
            y = np.array(conf.graphs_data.__getitem__(protocol_number=protocol.number, graph_type=self.graph_type).Y) * p

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

    def normalize(self, Y: np.ndarray, idx_min: Optional[int] = None, idx_max: Optional[int] = None) -> np.ndarray:
        if idx_min is not None and idx_max is not None:
            max_value = np.max(np.abs(Y[idx_min:idx_max]))
        else: 
            max_value = np.max(np.abs(Y))

        return Y / max_value

    
    def dump_XY_data(self, X: np.ndarray, Y: np.ndarray, fname: str) -> None:
        data = np.column_stack((X,Y))
        np.savetxt(fname, data, delimiter=' ')

    def check_conf(self, conf: Conformer, protocol: Protocol) -> bool:
        if not conf.active: 
            return False
        if not conf.graphs_data.__contains__(protocol.number):
            return False
        if not conf.graphs_data.__has_graph_type__(protocol.number, self.graph_type):
            return False
        return True 
    

    def diversity_function(self, a: np.ndarray, b: np.ndarray, w: Optional[np.ndarray] = None) -> float:
        # RMSD
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

        # after retrieving data, ensure we actually have peaks
        if self.energies.size == 0 or self.energies[self.energies!=0].size == 0:
            self.log.spectra_skip(self.graph_type)
            return

        if self.ref: 
            self.autoconvolution()
        else:
            self.SHIFT = self.defaults.shift
            self.FWHM = self.defaults.fwhm

            Y = self.convolute(energies=self.energies, impulses=self.impulse, shift=self.SHIFT, fwhm=self.FWHM)

            self.Y = self.normalize(Y)

            self.log.spectra_result(graph_type=self.graph_type, parameters={"Shift": self.SHIFT, "FWHM": self.FWHM}, msg=f"Using default parameters, Reference {self.graph_type} Spectra not found")

        if self.Y[~np.isnan(self.Y)].size > 0:
            self.log.debug(f'Saving {self.graph_type} spectra convoluted')
            self.dump_XY_data(self.X, self.Y, f'{self.graph_type}_p{self.protocol.number}_comp.xy')


    def autoconvolution(self) -> None:
        
        ref_norm = self.ref.Y

        def callback_optimizer(params: tuple) -> float:
            shift, fwhm = params
            st = datetime.now()
            Y_conv = self.convolute(self.energies, self.impulse, shift, fwhm)
            e1 = datetime.now()
            Y_conv = self.normalize(Y_conv, idx_min=self.ref.x_min_idx, idx_max=self.ref.x_max_idx)
            e2 = datetime.now()
            rmsd = self.diversity_function(Y_conv, ref_norm)
            e3 = datetime.now()
            # self.log.debug(f'{shift=:.2f}\t{fwhm=:.2f}\t{rmsd=:.2f}\t{e1-st}\t{e2-e1}\t{e2-st}\t{e3-st}')
            return rmsd
        
        initial_guess = [
            sum(self.shift_bounds)*.5, 
            sum(self.fwhm_bounds)*.5, 
        ]

        st = datetime.now()

        result = minimize(
            fun=callback_optimizer, x0=initial_guess, bounds=(self.shift_bounds, self.fwhm_bounds), options={"maxiter": 1000}#, method="Powell"
        ) 
        end = datetime.now()

        if result.success: 
            self.SHIFT, self.FWHM = result.x
            t = "Spectra convolution results:"
        else: 
            self.SHIFT, self.FWHM = self.defaults.shift, self.defaults.fwhm
            t = "Spectra convolution did NOT converged. Using default parameters:"


        Y = self.convolute(energies=self.energies, impulses=self.impulse, shift=self.SHIFT, fwhm=self.FWHM)
        self.Y = self.normalize(Y, idx_min=self.ref.x_min_idx, idx_max=self.ref.x_max_idx)

        diversity = self.diversity_function(self.Y[self.ref.x_min_idx:self.ref.x_max_idx], ref_norm[self.ref.x_min_idx:self.ref.x_max_idx])
        similarity = ((1 if self.graph_type not in CHIRALS else 2)-diversity)/(1 if self.graph_type not in CHIRALS else 2)*100

        diversity_unw = self.diversity_function(self.Y[self.ref.x_min_idx:self.ref.x_max_idx], ref_norm[self.ref.x_min_idx:self.ref.x_max_idx], w=np.ones_like(self.Y[self.ref.x_min_idx:self.ref.x_max_idx]))
        similarity_unw = ((1 if self.graph_type not in CHIRALS else 2)-diversity_unw)/(1 if self.graph_type not in CHIRALS else 2)*100

        self.log.spectra_result(graph_type=self.graph_type, parameters=
                                {"Shift": self.SHIFT, "FWHM": self.FWHM, "Similarity": similarity, "Similarity Unweighted": similarity_unw, "Time": (end-st), "Cycle": f"{result.nfev}"},
                                msg=t)



    def gaussian(self, x0: np.ndarray, I: np.ndarray, fwhm: float) -> np.ndarray:
        return gaussian_njit(self.X, x0, I, fwhm)
    
    def lorentzian(self, x0: np.ndarray, I: np.ndarray, fwhm: float) -> np.ndarray:
        return lorentzian_njit(self.X, x0, I, fwhm)





@njit(fastmath=True, cache=True)
def gaussian_njit(X, x0, I, fwhm):
    n_x = X.shape[0]
    n_peaks = x0.shape[0]
    Y = np.zeros(n_x)
    if n_peaks == 0:
        return Y

    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    norm = 1.0 / (sigma * np.sqrt(2 * np.pi))
    inv_sigma = 1.0 / sigma

    for j in prange(n_x):  # parallelizzato su X
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


@njit(fastmath=True, cache=True)
def lorentzian_njit(X, x0, I, fwhm):
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

@njit(fastmath=True, cache=True)
def diversity_function_njit(a, b, weight, max_val):
    diff = a - b
    s = 0.0
    n = diff.shape[0]
    for i in prange(n):
        s += weight[i] * diff[i] * diff[i]

    return np.sqrt(s / n) / max_val

