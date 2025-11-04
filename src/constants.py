from scipy.constants import R, c, h, electron_volt, Boltzmann, N_A
from scipy.constants import physical_constants


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


c = c * 100  # convert speed of light in cm/s
J_TO_H = physical_constants["joule-hartree relationship"][0]
AMU_TO_KG = physical_constants["atomic mass constant"][0]

EH_TO_KCAL = 627.5096080305927

CHIRALS = ["VCD", "ECD"]

GRAPHS = ['IR', 'VCD', 'UV', 'ECD']

CONVERT_B = {
    'GHz': 29.979000,
}