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
    """
    Format a frequency value as a string with 2 decimal places.

    Args:
        freq (float): The frequency value to format.

    Returns:
        str: Formatted string (e.g., "123.45").
    """

    return f"{freq:.2f}"


def get_conf_parameters(
    conf: Conformer, number: int, output: str, p, time: datetime, temp: float, log: Logger, linear: bool = False, cut_off: float = 100, alpha: float = 4, P: float = 101.325
) -> bool:
    """
    Extract and calculate parameters for a single conformer from the QM output.

    Parses energy, geometry, and frequencies. Calculates thermodynamic properties
    (G, H, S) using qRRHO if frequencies are available. Updates the conformer's
    EnergyStore and SpectralStore.

    Args:
        conf (Conformer): The conformer object to update.
        number (int): Protocol step number.
        output (str): Path to the output file.
        p (Protocol): Protocol configuration object.
        time (datetime): Timestamp/Duration of the calculation.
        temp (float): Temperature [K].
        log (Logger): Logger instance.
        linear (bool, optional): Treat molecule as linear. Defaults to False.
        cut_off (float, optional): qRRHO cut-off frequency [cm^-1]. Defaults to 100.
        alpha (float, optional): qRRHO damping factor. Defaults to 4.
        P (float, optional): Pressure [kPa]. Defaults to 101.325.

    Returns:
        bool: True if parsing was successful (even if calculation failed/diverged),
              False if a critical error occurred.

    Raises:
        IOError: If frequencies are expected but not found in the output.
    """

    try:
        parser = PARSER_REGISTRY[p.calculator](
            output_name=os.path.join(conf.folder, output), log=log)

        # if calculation crashed (detected by parser checks), skip conformer safely
        if not parser.correct_exiting:
            conf.active = False
            log.warning(f"Parser detected abnormal termination for Conf {conf.number}. Deactivating.")
            return True

        e = parser.parse_energy()

        if p.opt or 'opt' in p.add_input.lower():
            try:
                conf.last_geometry = parser.parse_geom().copy()
            except Exception as e:
                log.warning(f"Failed to parse geometry for Conf {conf.number}: {e}")

            if not parser.opt_done():
                if p.skip_opt_fail:
                    conf.active = False
                    log.warning(
                        f'{log.WARNING} Optimization did not correctly converge. Conf {conf.number} will be deactivated')
                    return True
                else: 
                    raise RuntimeError(f"Calculation for Conf {conf.number} did not finish: geometry not converged correctly")

        freq = np.array([])
        if p.freq or 'freq' in p.add_input.lower() or 'freq' in p.functional.lower():
            freq, ir, vcd = parser.parse_freq()
            if freq.size == 0:
                log.critical(
                    f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nNo frequency present in the calculation output for Conf {conf.number}.\n{'='*20}\n"
                )
                raise IOError("No frequency in the output file")

            freq *= p.freq_fact

        B_vec, M_vec = parser.parse_B_m()
        b = np.linalg.norm(B_vec) if B_vec is not None else 0
        m = np.linalg.norm(M_vec) if M_vec is not None else 0

        g = np.nan
        g_e, zpve, H, S = np.nan, np.nan, np.nan, np.nan
        
        if freq.size > 0:
            try:
                g, zpve, H, S = free_gibbs_energy(
                    SCF=e, T=temp, freq=freq, mw=conf.weight_mass, B=B_vec, m=p.mult,
                    linear = linear, cut_off = cut_off, alpha = alpha, P = P
                )
                H = H
                g_e = g - e
            except Exception as thermo_err:
                 log.error(f"Error computing RRHO thermodynamics for Conf {conf.number}: {thermo_err}")
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
                E=e if e else np.nan,
                G=g if not np.isnan(g) else np.nan,
                B=b if b else 1,
                m=m if m else 1,
                time=time,
                G_E=g_e if not np.isnan(g) and e else np.nan,
                zpve=zpve if not np.isnan(g) else np.nan,
                H=H if not np.isnan(g) else np.nan,
                S=S if not np.isnan(g) else np.nan,
                Freq=freq,
                B_vec=B_vec,
                m_vec=M_vec,
            )
        )
        log.debug(f'{log.TICK} Energy Data are stored correctly')

        try:
            _, ir, vcd = parser.parse_freq()
            uv, ecd = parser.parse_tddft()
            
            if ir.size == 0: ir = np.zeros((1, 2))
            if vcd.size == 0: vcd = np.zeros((1, 2))
            if uv.size == 0: uv = np.zeros((1, 2))
            if ecd.size == 0: ecd = np.zeros((1, 2))

            for label, graph in zip(GRAPHS, [ir, vcd, uv, ecd]):
                # Protezione su indici array vuoti
                if graph.shape[0] > 0 and graph.shape[1] >= 2:
                     conf.graphs_data.add(protocol_number=number, graph_type=label,
                                     record=SpectralRecord(X=graph[:, 0], Y=graph[:, 1]))
        except Exception as spec_err:
            log.warning(f"Error parsing spectral data for Conf {conf.number}: {spec_err}")

        log.debug(f'{log.TICK} Graphs Data are stored correctly')

        return True

    except Exception as e:
        # 3. CATTURA IL CRASH SILENZIOSO
        log.error(f"UNHANDLED EXCEPTION while parsing Conf {conf.number}: {str(e)}")
        import traceback
        log.debug(traceback.format_exc()) 
        
        conf.active = False 
        return False