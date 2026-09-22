import numpy as np
from typing import Optional
from ensemble_analyzer.constants import h, c, J_TO_H, Boltzmann, N_A

GRIMME_BAV = 1.00e-44  # kg * m^2 (default average moment of inertia, Grimme 2012)

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


def calc_zpe(frequency: Optional[np.ndarray] = None) -> float:
    r"""
    Calculate the Zero Point Energy.

    .. math::
        ZPE = \sum_{\nu}^\text{freq} \frac 12 h\nu c

    Args:
        frequency (np.ndarray, optional): Frequency list. Defaults to np.array([0]).

    Returns:
        float: Zero point energy in Eh.
    """
    if frequency is None or len(frequency) == 0:
        return 0.0
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
    freq: np.ndarray, T: float, cut_off: float = 100, alpha: int = 4
) -> float:
    r"""
    Harmonic vibrational thermal energy contribution (Grimme mRRHO style).

    .. math::
        U_{\text{vib}} = \sum_{\nu}^{\text{freq}} \frac{h\nu c}{e^{\frac{h\nu c}{k_bT}} - 1}

    In Grimme's 2012 mRRHO approach, the vibrational internal energy and enthalpy
    are computed with the standard harmonic oscillator partition function.

    Args:
        freq (np.ndarray): Vibrational frequencies [cm-1].
        T (float): Temperature [K].
        cut_off (float, optional): Frequency cutoff (retained for compatibility).
        alpha (int, optional): Damping factor (retained for compatibility).

    Returns:
        float: Vibrational thermal energy in Eh.
    """
    freq = np.asarray(freq, dtype=float)
    freq = freq[freq > 0]
    if len(freq) == 0:
        return 0.0
    return float(np.sum(calc_qRRHO_energy(freq, T)) * J_TO_H)


def calc_translational_entropy(MW: float, T: float, P: float) -> float:
    r"""
    Translational entropy.

    .. math::
        S_{trans} = k_b \left(\frac 52 + \ln\left(\sqrt{\frac{2πMWk_bT}{N_A*h^2}}^3 \frac {k_bT}{p}\right)\right)

    Args:
        MW (float): Molecular weight.
        T (float): Temperature.
        P (float): Pressure [kPa]. Defaults to 101.325.

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
        S_R &=& k_b \left(\ln\left(\frac{q_{rot}}{σ}\right) + 1.5\right)

    Args:
        B (np.array): Rotational constant [cm-1].
        T (float): Temperature.
        symno (int, optional): Number of symmetry, in relation of the Point Group of the molecule (σ). Defaults to 1.
        linear (bool, optional): If molecule is linear. Defaults to False.

    Returns:
        float: Rotational entropy in Eh.
    """
    B_arr = np.asarray(B, dtype=float)
    pos_B = B_arr[B_arr > 0]
    if len(pos_B) == 0:
        return 0.0

    rot_temperature = h * c * pos_B / Boltzmann

    if linear:
        qrot = T / rot_temperature[0]
    else:
        qrot = np.sqrt(np.pi * T**3 / np.prod(rot_temperature))

    symno = max(1, int(symno))
    return Boltzmann * (np.log(qrot / symno) + 1 + (0 if linear else 0.5)) * J_TO_H


def calc_S_V_grimme(freq: np.ndarray, T: float) -> np.ndarray:
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


def calc_S_R_grimme(freq: np.ndarray, T: float, B: np.ndarray | float | None = None) -> np.ndarray:
    r"""
    Free rotor entropy factor used for the quasi-RRHO interpolation (Grimme 2012).

    .. math::
        \mu &= \frac{h}{8\pi^2 \nu c}\\
        \mu' &= \frac{\mu B_{\text{av}}}{\mu + B_{\text{av}}}\\
        S_{\text{rot}} &= k_b \left( \frac{1}{2} + \ln\sqrt{\frac{8\pi^3 \mu' k_b T}{h^2}} \right)

    Args:
        freq (np.ndarray): Frequencies [cm-1].
        T (float): Temperature [K].
        B (np.ndarray, float, optional): Rotational constants [cm-1] or average moment of inertia [kg*m^2].
                                         If None or empty, GRIMME_BAV (1.00e-44 kg*m^2) is used.

    Returns:
        np.ndarray: Free rotor entropy in J/K.
    """
    freq = np.asarray(freq, dtype=float)
    if B is None:
        bav = GRIMME_BAV
    elif isinstance(B, (int, float, np.floating, np.integer)):
        bav = float(B) if B > 0 else GRIMME_BAV
    else:
        B_arr = np.asarray(B, dtype=float)
        pos_B = B_arr[B_arr > 0]
        if len(pos_B) > 0:
            # B_i in cm^-1; moment of inertia I_i = h / (8 * pi^2 * c * B_i) in kg * m^2
            moments = h / (8 * np.pi**2 * pos_B * c)
            bav = float(np.mean(moments))
        else:
            bav = GRIMME_BAV

    mu = h / (8 * np.pi**2 * freq * c)
    mu_prime = mu * bav / (mu + bav)
    f = 8 * np.pi**3 * mu_prime * Boltzmann * T / h**2

    return (0.5 + np.log(f**0.5)) * Boltzmann


def calc_vibrational_entropy(freq: np.ndarray, T: float, B: np.ndarray | float | None = None, cut_off=100, alpha=4) -> float:
    r"""
    Vibrational entropy using Grimme's quasi-RRHO interpolation model (2012).

    .. math::
        \sum_{\nu}^{freq} \left(w(\nu) S_{\text{vib,HO}}(\nu) + (1-w(\nu)) S_{\text{rot}}(\nu, T, B)\right)

    Args:
        freq (np.ndarray): Frequencies [cm-1].
        T (float): Temperature [K].
        B (np.ndarray, float, optional): Rotational constants [cm-1] or moment of inertia. Defaults to GRIMME_BAV.
        cut_off (float, optional): Cut off for the damping of the frequency. Defaults to 100.
        alpha (float, optional): Damping factor. Defaults to 4.

    Returns:
        float: Vibrational entropy [Eh/K].
    """
    freq = np.asarray(freq, dtype=float)
    freq = freq[freq > 0]
    if len(freq) == 0:
        return 0.0

    s_damp = calc_damp(freq, cut_off, alpha)
    return float(
        np.sum(
            calc_S_V_grimme(freq, T) * s_damp
            + (1 - s_damp) * calc_S_R_grimme(freq, T, B)
        )
        * J_TO_H
    )


def calc_vibrational_entropy_truhlar(
    freq: np.ndarray, T: float, cut_off: float = 100.0
) -> float:
    r"""
    Vibrational entropy calculated with Truhlar's quasi-harmonic cutoff model (2011).

    All frequencies below cut_off are raised to cut_off:
    .. math::
        \tilde{\omega}_i = \max(\omega_i, \omega_{\text{cut}})

    Args:
        freq (np.ndarray): Vibrational frequencies [cm-1].
        T (float): Temperature [K].
        cut_off (float, optional): Frequency cutoff in cm-1. Defaults to 100.0.

    Returns:
        float: Vibrational entropy in Eh/K.
    """
    freq = np.asarray(freq, dtype=float)
    freq = freq[freq > 0]
    if len(freq) == 0:
        return 0.0
    shifted_freq = np.maximum(freq, cut_off)
    return float(np.sum(calc_S_V_grimme(shifted_freq, T)) * J_TO_H)


def calc_electronic_entropy(m: int) -> float:
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
    cut_off: float = 100,
    alpha: int = 4,
    P: float = 101.325,
    symno: int = 1,
    model: str = "grimme",
) -> tuple[float, float, float, float]:
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
        symno (int, optional): Rotational symmetry number (sigma). Defaults to 1.
        model (str, optional): Quasi-RRHO model ('grimme' or 'truhlar'). Defaults to 'grimme'.

    Returns:
        tuple[float, float, float, float]: (G [Eh], zpve [Eh], h_corr [Eh], S [Eh/K])
    """
    freq = freq[freq > 0]

    zpve = calc_zpe(freq)

    U_trans = calc_translational_energy(T)
    U_rot = calc_rotational_energy(T, linear) if zpve > 0 else 0
    U_vib = calc_vibrational_energy(freq, T, cut_off, alpha)

    h_corr = zpve + U_trans + U_rot + U_vib + Boltzmann * T * J_TO_H
    H = SCF + h_corr

    S_elec = calc_electronic_entropy(m)
    if model.lower() == "truhlar":
        S_vib = calc_vibrational_entropy_truhlar(freq, T, cut_off=cut_off)
    elif model.lower() == "grimme":
        S_vib = calc_vibrational_entropy(freq, T, B, cut_off, alpha)
    else:
        raise ValueError(f"Unknown quasi-RRHO model '{model}'. Choose 'grimme' or 'truhlar'.")

    S_rot = calc_rotational_entropy(B, T, symno=symno, linear=linear)
    S_trans = calc_translational_entropy(mw, T, P)

    S = S_trans + S_rot + S_vib + S_elec

    return H - T * S, zpve, h_corr, S
