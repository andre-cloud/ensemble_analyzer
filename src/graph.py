import numpy as np
import os, sys
import pickle as pl

import matplotlib.pyplot as plt
from scipy.optimize import minimize
from typing import Optional, List, Union
import logging

logging.getLogger("matplotlib").disabled = False
logging.getLogger("numba").disabled = False

logging.getLogger("numba").setLevel(logging.WARNING)

try:
    from src.constants import *
    from src._spectral.experimental import ExperimentalGraph
    from src._spectral.comp_electronic import ComputedElectronic
    from src._spectral.comp_vibronic import ComputedVibronic
    from src._spectral.compare import ComparedGraph
    from src._protocol.protocol import Protocol
    from src._conformer.conformer import Conformer
except ModuleNotFoundError: 
    from constants import *
    from _spectral.experimental import ExperimentalGraph
    from _spectral.comp_electronic import ComputedElectronic
    from _spectral.comp_vibronic import ComputedVibronic
    from _spectral.compare import ComparedGraph
    from _protocol import Protocol
    from _conformer import Conformer


def eV_to_nm(eV):
    return FACTOR_EV_NM / eV

class_ = {
    'IR': ComputedVibronic, 
    'VCD': ComputedVibronic, 
    'UV': ComputedElectronic, 
    'ECD': ComputedElectronic, 
}




def main_spectra(ensemble: List[Conformer], protocol: Protocol, log: logging.Logger, invert: bool, interested_area: List[float], shift: Optional[Union[List[float],float,None]] = None, fwhm: Optional[Union[List[float],float,None]] = None, read_pop: Optional[str] = None, definition: Optional[int] = 4):
    
    for graph_type in list(class_.keys()):
        log.info("\n")
        ref = None
        fname = f"{graph_type.lower()}_ref.dat"
        if os.path.exists(os.path.join(os.getcwd(), fname)):
            ref = ExperimentalGraph(confs=[], protocol=protocol, graph_type=graph_type, log=log, interested_area=interested_area.get(VIBRO_OR_ELECTRO[graph_type]))
            ref.read()


        graph = class_[graph_type](
            confs=ensemble,
            graph_type=graph_type,
            log=log,
            protocol=protocol,
            ref=ref,
            invert=(graph_type in CHIRALS) and invert,
            shift_user=shift.get(VIBRO_OR_ELECTRO[graph_type], None),
            fwhm_user=fwhm.get(VIBRO_OR_ELECTRO[graph_type], None),
            read_population=read_pop,
            definition=definition
        )

        graph.compute_spectrum()

def plot_comparative_graphs(log, idxs=None, show=False, nm=True, show_ref_weight=False):
    for graph_type in GRAPHS:
        experimental_file = f"{graph_type.upper()}_ref_norm.xy" if os.path.exists(f"{graph_type.upper()}_ref_norm.xy") else None
        
        comp = ComparedGraph(graph_type=graph_type, experimental_file=experimental_file, log=log, protocol_index=idxs, nm=nm)
        if len(comp.data) > 0:
            comp.plot(show=show, show_ref_weight=show_ref_weight)