import numpy as np
try:
    from ensemble_analyzer.constants import *
except ModuleNotFoundError: 
    from ensemble_analyzer.constants import *

def calc_damp(frequency: np.ndarray, cut_off: float, alpha: int) -> np.ndarray:
    r"""
    Damping factor proportionate to frequency

    .. math::
        \frac {1}{1+(\frac {\text{cut_off}}{ν})^α}

    Damping factory NO measure unit

    :param frequency: Frequency list
    :type frequency: np.ndarray
    :param cut_off: cut off value, default is 100 cm-1
    :type cut_off: float
    :param alpha: dumping factor, default is 4
    :type alpha: int
    :return: dumping factor
    :rtype: np.ndarray
    """
    return 1 / (1 + (cut_off / frequency) ** alpha)


def calc_zpe(frequency: np.ndarray = np.array([0])) -> float:
    r"""Calculate the Zero Point Energy

    .. math::
        ZPE = \sum_{\nu}^\text{freq} \frac 12 h\nu c

    Return in Eh

    :param frequency: frequency list, defaults to np.array([0])
    :type frequency: np.ndarray, optional
    :return: zero point energy
    :rtype: float
    """
    return np.sum((h * frequency * c) / (2)) * J_TO_H


def calc_translational_energy(T: float) -> float:
    r"""Translational energy

    .. math::
        U_{trans} = \frac 32 K_bT

    Return in Eh

    :param T: temperature [K]
    :type T: float
    :return: translational energy
    :rtype: float
    """

    return 1.5 * Boltzmann * T * J_TO_H


def calc_rotational_energy(T: float, linear=False) -> float:
    r"""Rotational energy

    .. math::
        U_{rot} = \frac 32 K_bT\\
        U_{rot} = K_bT 

    Return in Eh

    :param T: temperature [K]
    :type T: float
    :param linear: if the molecule is linear, defaults to False
    :type linear: bool, optional
    :return: rotational energy
    :rtype: float
    """
    if linear:
        return Boltzmann * T * J_TO_H
    return 1.5 * Boltzmann * T * J_TO_H


def calc_qRRHO_energy(freq: np.ndarray, T: float) -> np.ndarray:
    r"""quasi-Rigid Rotor Harmonic Oscillator energy

    .. math::
        U = h\nu c \frac { e^{-\frac {h\nu c}{k_bT}} }{1-e^{-\frac {h\nu c}{k_bT}}}


    :param freq: frequency list
    :type freq: np.ndarray
    :param T: temperature [K]
    :type T: float
    :return: vibrational energy for each vibrational mode in Joule
    :rtype: np.ndarray
    """
    f = h * freq * c / (Boltzmann * T)
    return h * freq * c * np.exp(-f) / (1 - np.exp(-f))


def calc_vibrational_energy(
    freq: np.ndarray, T: float, cut_off: float, alpha: int
) -> float:
    r"""
    Vibrational energy calculated with qRRHO.

    .. math::
        \sum_{\nu}^{freq} \left( d H_{qRRHO}(freq, T) + (1 - d)k_bT\frac 12 \right)

    :param freq: frequency
    :type freq: array
    :param T: temperature
    :type T: float
    :param alpha: damping factor, default and unchangeable value is 4.
    :type alpha: int, optional
    :param cut_off: damping frequency, default 100 cm-1
    :type cut_off: float, optional

    :return: vibrational energy in Eh
    :rtype: float
    """
    h_damp = calc_damp(freq, cut_off=cut_off, alpha=alpha)
    return (
        np.sum(h_damp * calc_qRRHO_energy(freq, T) + (1 - h_damp) * Boltzmann * T * 0.5)
        * J_TO_H
    )


def calc_translational_entropy(MW: float, T: float, P: float) -> float:
    r"""
    Translational entropy

    .. math::
        S_{trans} = k_b \left(\frac 52 + \ln\left(\sqrt{\frac{2πMWk_bT}{N_A*h^2}}^3 \frac {k_bT}{p}\right)\right)


    :param MW: molecular weight
    :type MW: float
    :param T: temperature
    :type T: float
    :param P: pressure [Pa]
    :type P: float
    :param solv: solvent string
    :type solv: str


    :return: translational entropy in Eh
    :rtype: float
    """

    lambda_ = np.sqrt((2 * np.pi * MW * Boltzmann * T) / (1000 * N_A * h**2))
    V = (Boltzmann * T) / (P * 1000)

    return Boltzmann * (5 / 2 + np.log(lambda_**3 * V)) * J_TO_H


def calc_rotational_entropy(B, T, symno: int = 1, linear: bool = False) -> float:
    r"""
    Rotational entropy

    .. math::
        θ_R &=& \frac {hcB}{k_b}\\
        q_{rot} &=& \sqrt{\frac {πT^3}{θ_{Rx}θ_{Ry}θ_{Rz}}}\\
        S_R &=& k_b \left(\frac {\ln(q_{rot}}{σ} + 1.5\right)

    :param B: rotational constant [cm-1]
    :type B: np.array
    :param symno: number of symmetry, in relation of the Point Group of the molecule (σ), default is 1
    :type symno: int
    :param linear: if molecule is linear, default is False
    :type linear: bool

    :return: rotational entropy in Eh
    :rtype: float
    """
    rot_temperature = h * c * B / Boltzmann

    if linear:
        qrot = T / rot_temperature[0]
    else:
        qrot = np.sqrt(np.pi * T**3 / np.prod(rot_temperature))

    return Boltzmann * (np.log(qrot / symno) + 1 + (0 if linear else 0.5)) * J_TO_H


def calc_S_V_grimme(freq: np.array, T) -> np.array:
    r"""
    V factor used for the damping of the frequency

    .. math::
        V = \frac {\frac {hc\nu}{k_bT} k_b}{e^{\frac {hc\nu}{k_bT}} - 1} - k_b \ln\left(1 - e^{-\frac {hc\nu}{k_bT}}\right)

    :param freq: frequencies [cm-1]
    :type freq: np.array
    :param T: temperature [K]
    :type T: float

    :return: V factor in J
    :rtype: float
    """
    f = h * freq * c / (Boltzmann * T)
    return (f * Boltzmann) / (np.exp(f) - 1) - Boltzmann * np.log(1 - np.exp(-f))


def calc_S_R_grimme(freq: np.array, T: float, B: np.array) -> np.array:
    r"""
    R factor used for the damping of the frequency

    .. math::
        R = \frac 12 \left( 1+ \ln\left( \frac {8π^3 \frac {h}{8π^2\nu c} B k_bT} {\left(\frac {h}{8π^2\nu c}+B\right)h^2} \right)\right) k_b


    :param freq: frequencies [cm-1]
    :type freq: np.array
    :param T: temperature [K]
    :type T: float
    :param B: rotatory constant [cm-1]
    :type B: np.array

    :return: R factor in J
    :rtype: float
    """

    B = (np.sum(B * c) / len(B)) ** -1 * h
    mu = h / (8 * np.pi**2 * freq * c)
    f = 8 * np.pi**3 * (mu * B / (mu + B)) * Boltzmann * T / h**2

    return (0.5 + np.log(f**0.5)) * Boltzmann


def calc_vibrational_entropy(freq, T, B, cut_off=100, alpha=4) -> float:
    r"""
    Vibrational entropy

    .. math::
        \sum_{\nu}^{freq} \left(dV(\nu) + (1-d)R(\nu, T, B)\right)

    in formula :math:`d` is the dumping function

    :param freq: frequencies [cm-1]
    :type freq: list
    :param T: temperature [K]
    :type T: float
    :param B: rotational constant [cm-1]
    :type B: np.array
    :param cut_off: cut off for the damping of the frequency
    :type cut_off: float
    :param alpha: damping factor
    :type alpha: float

    :return: vibrational entropy [Eh]
    :rtype: float
    """

    s_damp = calc_damp(freq, cut_off, alpha)
    return (
        np.sum(
            calc_S_V_grimme(freq, T) * s_damp
            + (1 - s_damp) * calc_S_R_grimme(freq, T, B)
        )
        * J_TO_H
    )


def calc_electronic_entropy(m) -> float:
    r"""
    Electronic entropy

    .. math::
        S_{el} = k_b \ln(m)

    :param m: electronic multiplicity
    :type m: int

    :return: electronic entropy in Eh
    :rtype: float
    """
    return Boltzmann * np.log(m) * J_TO_H


def free_gibbs_energy(
    SCF: float,
    T: float,
    freq: np.ndarray,
    mw: float,
    B: np.ndarray,
    m: int,
    # defaults
    linear: bool = False,
    cut_off=100,
    alpha=4,
    P: float = 101.325,
) -> float:
    r"""
    Calculate Gibbs energy

    .. math::
        H &=& SCF + ZPVE + U_{trans} + U_{rot} + U_{vib} + k_bT\\
        S &=& S_{trans} + S_{rot} + S_{vib} + S_{el}\\
        G &=& H - TS

    :param SCF: self consistent field energy [Eh] + dispersions
    :type SCF: float
    :param T: temperature [K]
    :type T: float
    :param P: pressure [kPa], default is 101.325
    :type P: float
    :param B: rotational constant [cm-1]
    :type B: np.array
    :param m: spin multiplicity
    :type m: int

    :param linear: if molecule is linear
    :type linear: bool
    :param cut_off: frequency cut_off
    :type cut_off: float
    :param alpha: frequency damping factor
    :type alpha: int
    :param P: pressure [kPa]
    :type P: float

    :return: Gibbs energy
    :rtype: float
    """
    freq = freq[freq > 0]

    zpve = calc_zpe(freq)

    U_trans = calc_translational_energy(T)
    U_rot = calc_rotational_energy(T, linear) if zpve > 0 else 0
    U_vib = calc_vibrational_energy(freq, T, cut_off, alpha)

    h = zpve + U_trans + U_rot + U_vib + Boltzmann * T * J_TO_H
    H = SCF + h

    S_elec = calc_electronic_entropy(m)
    S_vib = calc_vibrational_entropy(freq, T, B, cut_off, alpha)
    S_rot = calc_rotational_entropy(B, T, linear=linear)
    S_trans = calc_translational_entropy(mw, T, P)

    S = S_trans + S_rot + S_vib + S_elec

    return H - T * S, zpve, h, S


if __name__ == "__main__":
    from ensemble_analyzer.parser_parameter import get_param, get_freq
    import sys

    args = sys.argv[1:]
    *output, calc, T = args

    for i in output:
        with open(i) as f:
            fl = f.readlines()

        e = float(
            list(filter(lambda x: get_param(x, calc, "E"), fl))[-1].strip().split()[-1]
        )

        B = np.array(
            list(filter(lambda x: get_param(x, calc, "B"), fl))[-1]
            .strip()
            .split(":")[-1]
            .split(),
            dtype=float,
        )

        mw = float(
            [i for i in fl if "Total Mass" in i][0]
            .strip()
            .split("...")[1]
            .split()[0]
            .strip()
        )

        freq = get_freq(fl, calc)
        im_freq = freq[freq < 0]

        g = free_gibbs_energy(SCF=e, T=float(T), freq=freq[freq > 0], mw=mw, B=B, m=1)
        print(
            f'{i} --- G with mRRHO @ T={T}: {g} Eh     Calculation ended with {len(im_freq)} imaginary frequencies {" ".join(list(map(str, im_freq))) if len(im_freq) > 0 else ""}'
        )
