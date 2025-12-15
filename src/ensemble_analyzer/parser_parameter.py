import re
import os
import numpy as np

from ensemble_analyzer.constants import regex_parsing
from ensemble_analyzer.rrho import free_gibbs_energy
from ensemble_analyzer.constants import *
from ensemble_analyzer._parsers.base import PARSER_REGISTRY

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._conformer.energy_data import EnergyRecord
from ensemble_analyzer._conformer.spectral_data import SpectralRecord

from ensemble_analyzer._logger.logger import Logger

from datetime import datetime


def tranform_float(freq) -> str:
    """Transform into a float number a string. Thought for a map function

    Args:
        freq (float): Frequency
    """
    return f"{freq:.2f}"


def get_conf_parameters(
    conf: Conformer, number: int, output: str, p, time: datetime, temp: float, log: Logger, linear: bool = False, cut_off: float = 100, alpha: float = 4, P: float = 101.325
) -> bool:
    """Obtain the parameters for a conformer: E, G, B, m

    Args:
        conf (Conformer): Conformer
        number (int): Protocol number
        output (str): Output File
        p (Protocol): Protocol
        time (datetime): Time of execution
        temp (float): Temperature [K]
        log (Logger): Logging system

    Raises:
        IOError: When frequency are not found

    Returns:
        bool: Correct parsing
    """

    parser = PARSER_REGISTRY[p.calculator](
        output_name=os.path.join(conf.folder, output), log=log)

    # if calculation crashed, skip conformer
    if not parser.correct_exiting:
        conf.active = False
        return True

    e = parser.parse_energy()

    if p.opt or 'opt' in p.add_input.lower():
        conf.last_geometry = parser.parse_geom().copy()

        if not parser.opt_done():
            if p.skip_opt_fail:
                conf.active = False
                log.warning(
                    f'{log.WARNING} Optimization did not correctly converge (maybe increase number of iteration). Conf {conf.number} will be deactivated')
                return True
            else: 
                raise f"Calculation for Conf {conf.number} did not finished: geometry not converged correctly"

        # TODO: LOGICA PER UN'OTTIMIZZAZIONE NON COMPLETATA
        # Si potrebbe rilanciare l'ottimizzazione...

    freq = np.array([])
    if p.freq or 'freq' in p.add_input.lower() or 'freq' in p.functional.lower():
        freq, ir, vcd = parser.parse_freq()
        if freq.size == 0:
            # TODO: creare un logger error per questa cosa
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nNo frequency present in the calculation output.\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise IOError("No frequency in the output file")

        freq *= p.freq_fact

    B_vec, M_vec = parser.parse_B_m()
    b = np.linalg.norm(B_vec)
    m = np.linalg.norm(M_vec)

    g = np.nan
    g_e, zpve, H, S = np.nan, np.nan, np.nan, np.nan
    if freq.size > 0:
        g, zpve, H, S = free_gibbs_energy(
            SCF=e, T=temp, freq=freq, mw=conf.weight_mass, B=B_vec, m=p.mult,
            linear = linear, cut_off = cut_off, alpha = alpha, P = P
        )
        H = H
        g_e = g - e
    else:
        prev_energies = conf.energies.__getitem__(int(number) - 1)
        g_e = prev_energies.G_E

        if not np.isnan(g_e):
            g = e + g_e
            zpve = prev_energies.zpve
            H = prev_energies.H
            S = prev_energies.S
        else:
            log.missing_previous_thermo(conformer_id=conf.number)

    conf.energies.add(
        number,
        EnergyRecord(
            E=e if e else e,  # Electronic Energy [Eh]
            G=g if not np.isnan(g) else np.nan,  # Free Gibbs Energy [Eh]
            B=b if b else 1,  # Rotatory Constant [cm-1]
            m=m if m else 1,  # dipole momenti [Debye]
            time=time,  # elapsed time [sec]
            G_E=g_e if not np.isnan(g) and e else np.nan,  # G-E [Eh]
            zpve=zpve if not np.isnan(g) else np.nan,  # Zero Point Energy [Eh]
            H=H if not np.isnan(g) else np.nan,  # Enthalpy correction [Eh]
            S=S if not np.isnan(g) else np.nan,  # Entropy [Eh],
            Freq=freq,  # Frequencies
            B_vec=B_vec,  # rotational vector constant [cm-1]
            m_vec=M_vec,  # Dipole moment vector [Debye]
        )
    )
    log.debug(f'{log.TICK} Energy Data are stored correctly')

    _, ir, vcd = parser.parse_freq()
    uv, ecd = parser.parse_tddft()
    for label, graph in zip(GRAPHS, [ir, vcd, uv, ecd]):
        conf.graphs_data.add(protocol_number=number, graph_type=label,
                             record=SpectralRecord(X=graph[:, 0], Y=graph[:, 1]))
    log.debug(f'{log.TICK} Graphs Data are stored correctly')

    return True
