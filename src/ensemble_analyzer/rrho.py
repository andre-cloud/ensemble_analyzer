import numpy as np
try:
    from ensemble_analyzer.constants import *
except ModuleNotFoundError: 
    from ensemble_analyzer.constants import *

def calc_damp(frequency: np.ndarray, cut_off: float, alpha: int) -> np.ndarray:
    r"""
    Damping factor proportionate to frequency.

    .. math::
        \frac {1}{1+(\frac {\text{cut_off}}{ν})^α}

    Damping factor has NO measure unit.

    Args:
        frequency (np.ndarray): Frequency list.
        cut_off (float): Cut off value, default is 100 cm-1.
        alpha (int): Damping factor, default is 4.

    Returns:
        np.ndarray: Damping factor.
    """
    return 1 / (1 + (cut_off / frequency) ** alpha)


def calc_zpe(frequency: np.ndarray = np.array([0])) -> float:
    r"""
    Calculate the Zero Point Energy.

    .. math::
        ZPE = \sum_{\nu}^\text{freq} \frac 12 h\nu c

    Args:
        frequency (np.ndarray, optional): Frequency list. Defaults to np.array([0]).

    Returns:
        float: Zero point energy in Eh.
    """
    return np.sum((h * frequency * c) / (2)) * J_TO_H


def calc_translational_energy(T: float) -> float:
    r"""
    Translational energy.

    .. math::
        U_{trans} = \frac 32 K_bT

    Args:
        T (float): Temperature [K].

    Returns:
        float: Translational energy in Eh.
    """
    return 1.5 * Boltzmann * T * J_TO_H


def calc_rotational_energy(T: float, linear=False) -> float:
    r"""
    Rotational energy.

    .. math::
        U_{rot} = \frac 32 K_bT\\
        U_{rot} = K_bT 

    Args:
        T (float): Temperature [K].
        linear (bool, optional): If the molecule is linear. Defaults to False.

    Returns:
        float: Rotational energy in Eh.
    """
    if linear:
        return Boltzmann * T * J_TO_H
    return 1.5 * Boltzmann * T * J_TO_H


def calc_qRRHO_energy(freq: np.ndarray, T: float) -> np.ndarray:
    r"""
    quasi-Rigid Rotor Harmonic Oscillator energy.

    .. math::
        U = h\nu c \frac { e^{-\frac {h\nu c}{k_bT}} }{1-e^{-\frac {h\nu c}{k_bT}}}

    Args:
        freq (np.ndarray): Frequency list.
        T (float): Temperature [K].

    Returns:
        np.ndarray: Vibrational energy for each vibrational mode in Joule.
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

    Args:
        freq (np.ndarray): Frequency array.
        T (float): Temperature.
        cut_off (float): Damping frequency, default 100 cm-1.
        alpha (int): Damping factor, default and unchangeable value is 4.

    Returns:
        float: Vibrational energy in Eh.
    """
    h_damp = calc_damp(freq, cut_off=cut_off, alpha=alpha)
    return (
        np.sum(h_damp * calc_qRRHO_energy(freq, T) + (1 - h_damp) * Boltzmann * T * 0.5)
        * J_TO_H
    )


def calc_translational_entropy(MW: float, T: float, P: float) -> float:
    r"""
    Translational entropy.

    .. math::
        S_{trans} = k_b \left(\frac 52 + \ln\left(\sqrt{\frac{2πMWk_bT}{N_A*h^2}}^3 \frac {k_bT}{p}\right)\right)

    Args:
        MW (float): Molecular weight.
        T (float): Temperature.
        P (float): Pressure [Pa].

    Returns:
        float: Translational entropy in Eh.
    """

    lambda_ = np.sqrt((2 * np.pi * MW * Boltzmann * T) / (1000 * N_A * h**2))
    V = (Boltzmann * T) / (P * 1000)

    return Boltzmann * (5 / 2 + np.log(lambda_**3 * V)) * J_TO_H


def calc_rotational_entropy(B, T, symno: int = 1, linear: bool = False) -> float:
    r"""
    Rotational entropy.

    .. math::
        θ_R &=& \frac {hcB}{k_b}\\
        q_{rot} &=& \sqrt{\frac {πT^3}{θ_{Rx}θ_{Ry}θ_{Rz}}}\\
        S_R &=& k_b \left(\frac {\ln(q_{rot}}{σ} + 1.5\right)

    Args:
        B (np.array): Rotational constant [cm-1].
        T (float): Temperature.
        symno (int, optional): Number of symmetry, in relation of the Point Group of the molecule (σ). Defaults to 1.
        linear (bool, optional): If molecule is linear. Defaults to False.

    Returns:
        float: Rotational entropy in Eh.
    """
    rot_temperature = h * c * B / Boltzmann

    if linear:
        qrot = T / rot_temperature[0]
    else:
        qrot = np.sqrt(np.pi * T**3 / np.prod(rot_temperature))

    return Boltzmann * (np.log(qrot / symno) + 1 + (0 if linear else 0.5)) * J_TO_H


def calc_S_V_grimme(freq: np.array, T) -> np.array:
    r"""
    V factor used for the damping of the frequency.

    .. math::
        V = \frac {\frac {hc\nu}{k_bT} k_b}{e^{\frac {hc\nu}{k_bT}} - 1} - k_b \ln\left(1 - e^{-\frac {hc\nu}{k_bT}}\right)

    Args:
        freq (np.array): Frequencies [cm-1].
        T (float): Temperature [K].

    Returns:
        np.array: V factor in J.
    """
    f = h * freq * c / (Boltzmann * T)
    return (f * Boltzmann) / (np.exp(f) - 1) - Boltzmann * np.log(1 - np.exp(-f))


def calc_S_R_grimme(freq: np.array, T: float, B: np.array) -> np.array:
    r"""
    R factor used for the damping of the frequency.

    .. math::
        R = \frac 12 \left( 1+ \ln\left( \frac {8π^3 \frac {h}{8π^2\nu c} B k_bT} {\left(\frac {h}{8π^2\nu c}+B\right)h^2} \right)\right) k_b

    Args:
        freq (np.array): Frequencies [cm-1].
        T (float): Temperature [K].
        B (np.array): Rotatory constant [cm-1].

    Returns:
        np.array: R factor in J.
    """

    B = (np.sum(B * c) / len(B)) ** -1 * h
    mu = h / (8 * np.pi**2 * freq * c)
    f = 8 * np.pi**3 * (mu * B / (mu + B)) * Boltzmann * T / h**2

    return (0.5 + np.log(f**0.5)) * Boltzmann


def calc_vibrational_entropy(freq, T, B, cut_off=100, alpha=4) -> float:
    r"""
    Vibrational entropy.

    .. math::
        \sum_{\nu}^{freq} \left(dV(\nu) + (1-d)R(\nu, T, B)\right)

    In formula :math:`d` is the dumping function.

    Args:
        freq (list): Frequencies [cm-1].
        T (float): Temperature [K].
        B (np.array): Rotational constant [cm-1].
        cut_off (float, optional): Cut off for the damping of the frequency. Defaults to 100.
        alpha (float, optional): Damping factor. Defaults to 4.

    Returns:
        float: Vibrational entropy [Eh].
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
    Electronic entropy.

    .. math::
        S_{el} = k_b \ln(m)

    Args:
        m (int): Electronic multiplicity.

    Returns:
        float: Electronic entropy in Eh.
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
    Calculate Gibbs energy.

    .. math::
        H &=& SCF + ZPVE + U_{trans} + U_{rot} + U_{vib} + k_bT\\
        S &=& S_{trans} + S_{rot} + S_{vib} + S_{el}\\
        G &=& H - TS

    Args:
        SCF (float): Self consistent field energy [Eh] + dispersions.
        T (float): Temperature [K].
        freq (np.ndarray): Frequencies array.
        mw (float): Molecular weight.
        B (np.array): Rotational constant [cm-1].
        m (int): Spin multiplicity.
        linear (bool, optional): If molecule is linear. Defaults to False.
        cut_off (float, optional): Frequency cut_off. Defaults to 100.
        alpha (int, optional): Frequency damping factor. Defaults to 4.
        P (float, optional): Pressure [kPa]. Defaults to 101.325.

    Returns:
        float: Gibbs energy.
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
