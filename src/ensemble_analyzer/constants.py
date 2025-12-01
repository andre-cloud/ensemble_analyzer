from scipy.constants import R, c, h, electron_volt, Boltzmann, N_A
from scipy.constants import physical_constants
import numpy as np

import os

DEBUG = bool(os.getenv("DEBUG"))
MAX_TRY = 5


# Physical CONSTANTS

# R = 8.314462618 J/(mol K)
# h = 6.62607015e-34 J*s
# c = 2.9979245800E+10 cm/s
# Boltzmann = 1.380649e-23 J/K
# J_TO_H = 2.2937122783963e+17 Eh/J
# AMU_TO_KG = 1.6605390666e-27 kg*mol/g

FACTOR_EV_NM = h * c / (10**-9 * electron_volt)
FACTOR_EV_CM_1 = 1 / 8065.544  # to yield eV


def eV_to_nm(eV):
    eV = np.maximum(eV, 1e-2)
    return FACTOR_EV_NM / eV


c = c * 100  # convert speed of light in cm/s
J_TO_H = physical_constants["joule-hartree relationship"][0]
AMU_TO_KG = physical_constants["atomic mass constant"][0]

EH_TO_KCAL = 627.5096080305927
CAL_TO_J = 4.186


CHIRALS = ["VCD", "ECD"]

GRAPHS = ['IR', 'VCD', 'UV', 'ECD']

CONVERT_B = {
    'GHz': 29.979000,
}

VIBRO_OR_ELECTRO = {
    'IR': 'vibro',
    'VCD': 'vibro',
    'UV': 'electro',
    'ECD': 'electro',
}


MIN_RETENTION_RATE = 0.3


# Logger constants
LOG_FORMAT = "%(message)s"


def ordinal(n):
    return "%d-%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10:: 4])


regex_parsing = {
    "orca": {"ext": "out", },
    "gaussian": {"ext": "log", },
}


MARKERS = [
    ".", ",", "o", "v", "^", "<", ">", "1", "2", "3",
    "4", "8", "s", "p", "*", "h", "H", "+", "x", "D",
    "d", "|", "_", "P", "X",
]

DEFAULT_RESOLUTION = 500        # Grid resolution for contour plots
MIN_CONFORMERS_FOR_PCA = 50     # Minimum conformers needed for meaningful PCA
MIN_WEIGHTED_VALUE = 0.15       # Minimum value for the weight of the autoconvolution