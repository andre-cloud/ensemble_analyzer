from dataclasses import dataclass
import os, sys
import numpy as np

from typing import Union, Optional, Literal

from src._spectral.graph_default import GraphDefault
from src._spectral.base import BaseGraph

from src.constants import *

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