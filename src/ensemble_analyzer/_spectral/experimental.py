from dataclasses import dataclass
import os, sys
import numpy as np

from typing import Union, Optional, Literal

from ensemble_analyzer._spectral.graph_default import GraphDefault
from ensemble_analyzer._spectral.base import BaseGraph

from ensemble_analyzer.constants import *

@dataclass
class ExperimentalGraph(BaseGraph): 

    graph_type: Literal['IR', 'VCD', 'UV', 'ECD']

    def __post_init__(self):
        self.defaults = GraphDefault(self.graph_type)

        X = np.linspace(self.defaults.start, self.defaults.end, num=10**self.definition)
        self.X = X[np.argsort(X)]

        self.fname: Optional[str] = self.defaults.experimental_fname

    def read(self):

        self.load_file_experimental()
        self.process()
        self.y = self.normalize(self.y)
        self.Y = self.interpolate()
        self.calc_weighting_function()

        self.dump_XY_data(self.X, self.weight, f'{self.graph_type}_weighted.xy')

        self.dump_XY_data(self.X, self.Y, f'{self.graph_type}_ref_norm.xy')
        np.savetxt(f'{self.graph_type.upper()}_index_lim.xy', np.array([self.x_min_idx, self.x_max_idx]))

    def load_file_experimental(self):

        fname = os.path.join(os.getcwd(), self.fname)
        self.log.debug("Reading the reference data from " + fname)

        data = np.loadtxt(fname, dtype=np.float64)
        data = data[np.argsort(data[:, 0])]
        X, Y = np.hsplit(data, 2)

        self.x = X.ravel()
        self.y = self.normalize(Y.ravel())

    def interpolate(self):
        Y =  np.interp(self.X, self.x, self.y, left=0, right=0)
        return Y


    def process(self): 

        convert = False
        if self.graph_type in ['UV', 'ECD'] and np.min(self.x) > 20: 
            self.x = FACTOR_EV_NM / self.x
            convert = True
            if convert:
                self.x = self.x[::-1]
                self.y = self.y[::-1]
        
        self.x_min = float(np.min(self.x))
        self.x_max = float(np.max(self.x))

        self.x_min_idx = int(np.argmin((self.X - self.x_min)<0))
        self.x_max_idx = int(np.argmax((self.X - self.x_max)>0))
        # if convert: 
        #     self.x_min_idx, self.x_max_idx = self.x_max_idx, self.x_min_idx

    def calc_weighting_function(self): 

        LIM_EXP1, LIM_EXP2 = self.X[self.x_min_idx], self.X[self.x_max_idx]
        self.weight = np.zeros_like(self.X)
        ia = self.interested_area


        if ia is None: 
            self.weight[(self.X >= LIM_EXP1) & (self.X <= LIM_EXP2)] = 1
            return
        
        if ia is True: 
            INT1, INT2 = self.defaults.interested_area

        elif isinstance(ia, (list, tuple)) and len(ia) == 2:
            INT1, INT2 = sorted(self.interested_area)

        elif isinstance(ia, (int, float)): 
            INT1, INT2 = self.interested_area, self.interested_area
        else: 
            raise ValueError("Invalid format for interested_area")

        sigma1 = np.abs(LIM_EXP1-INT1) * .1
        sigma2 = np.abs(LIM_EXP2-INT2) * .1

        gau1 = self.gau(INT1, sigma1)
        gau1 /= np.max(gau1)
        gau2 = self.gau(INT2, sigma2)
        gau2 /= np.max(gau2)

        self.weight[(LIM_EXP1 < self.X) & (self.X < INT1)] = gau1[(LIM_EXP1 < self.X) & (self.X < INT1)]
        self.weight[(INT1 <= self.X) & (self.X <= INT2)] = 1
        self.weight[(INT2 < self.X) & (self.X < LIM_EXP2)] = gau2[(INT2 < self.X) & (self.X < LIM_EXP2)]
        self.weight[self.weight < MIN_WEIGHTED_VALUE] = MIN_WEIGHTED_VALUE

        self.weight[self.X < LIM_EXP1] = 0
        self.weight[self.X > LIM_EXP2] = 0

        return


    def _calc_sigma(self, x1: float, x0: float, v:Optional[float]=0.2, A:Optional[float]=1.):
        # x1 Reference Limit
        # x0 Interest boundary Limit
        return abs(x1 - x0) / np.sqrt(2 * np.log(A/v))
    
    def gau(self, x0, sigma):
        return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((self.X - x0)**2)/sigma**2)


