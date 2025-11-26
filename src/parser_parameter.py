import re
import os
import numpy as np

try:
    from src.regex_parsing import regex_parsing
    from src.rrho import free_gibbs_energy
    from src.constants import * 
    from src._parsers.base import PARSER_REGISTRY

except ImportError as e:  # pragma: no cover
    print(e)
    from regex_parsing import regex_parsing
    from rrho import free_gibbs_energy
    from constants import * 
    from _parsers.base import PARSER_REGISTRY


def tranform_float(freq):
    """Transform into a float number a string. Thought for a map function

    :param freq: frequency to be transformed
    :type freq: float

    :return: frequency transformed
    :rtype: str
    """
    return f"{freq:.2f}"


def get_conf_parameters(
    conf, number: int, output: str, p, time, temp: float, log
) -> bool:
    """
    Obtain the parameters for a conformer: E, G, B, m

    :param conf: conformer
    :type conf: Conformer
    :param number: protocol number
    :type number: int
    :param p: protocol executed
    :type p: Protocol
    :param time: elapsed time requested for the calculation
    :type time: datetime
    :param temp: temperature [K]
    :type temp: float
    :param log: logger instance
    :type log: logger

    :return: calculation ended correctly and not crashed due to server error
    :rtype: bool
    """

    parser = PARSER_REGISTRY[p.calculator](output_name=os.path.join(conf.folder, output), log=log)

    # if calculation crashed, skip conformer
    if not parser.correct_exiting:
        conf.active = False
        return True
    
    e = parser.parse_energy()
    
    if p.opt or 'opt' in p.add_input.lower() or 'opt' in p.functional.lower():
        conf.last_geometry = parser.parse_geom().copy()
        
        if p.skip_opt_fail: 
            if not parser.opt_done():
                conf.active = False

                return True
            
        # TODO: LOGICA PER UN'OTTIMIZZAZIONE NON COMPLETATA
        # Si potrebbe rilanciare l'ottimizzazione... 

    freq= np.array([])
    if p.freq or 'freq' in p.add_input.lower() or 'freq' in p.functional.lower():
        freq, ir, vcd = parser.parse_freq()
        if freq.size == 0: 
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nNo frequency present in the calculation output.\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise IOError("No frequency in the output file")
         
        freq *= p.freq_fact
        log.info(
            f"\tConf {conf.number} has {freq[freq < 0].size} imaginary frequency(s): {', '.join(list(map(tranform_float, freq[freq < 0])))}"
        )
        if freq[freq < 0].size > 0:
            log.info(
                f"\tThese are excluded from qRRHO calculation."
            )

    B_vec, M_vec = parser.parse_B_m()
    b = np.linalg.norm(B_vec)
    m = np.linalg.norm(M_vec)

    g = np.nan
    g_e, zpve, H, S = np.nan, np.nan, np.nan, np.nan
    if freq.size > 0:
        g, zpve, H, S = free_gibbs_energy(
            SCF=e, T=temp, freq=freq, mw=conf.weight_mass, B=B_vec, m=p.mult
        )
        H = H
        g_e = g - e
    else:
        prev_energies = conf.energies.get(str(int(number) - 1), {})
        g_e = prev_energies.get('G-E',np.nan)

        if not np.isnan(g_e):
            g = e + g_e
            zpve = prev_energies.get("zpve", np.nan)
            H = prev_energies.get("H", np.nan)
            S = prev_energies.get("S", np.nan)
        else:
            log.missing_previous_thermo(conformer_id = conf.number)

    conf.energies[str(number)] = {
        "E": e if e else e,  # Electronic Energy [Eh]
        "G": g if not np.isnan(g) else np.nan,  # Free Gibbs Energy [Eh]
        "B": b if b else 1,  # Rotatory Constant [cm-1]
        "m": m if m else 1,  # dipole momenti [Debye]
        "time": time,  # elapsed time [sec]
        "G-E": g_e if not np.isnan(g) and e else np.nan,  # G-E [Eh]
        "zpve": zpve if not np.isnan(g) else np.nan,  # Zero Point Energy [Eh]
        "H": H if not np.isnan(g) else np.nan,  # Enthalpy correction [Eh]
        "S": S if not np.isnan(g) else np.nan,  # Entropy [Eh]
    }

    freq, ir, vcd = parser.parse_freq()
    uv, ecd = parser.parse_tddft()
    conf.energies[str(number)]["graph"] = {}
    for label, graph in zip(GRAPHS, [ir, vcd, uv, ecd]):
        conf.energies[str(number)]["graph"][label] = {'x': graph[:,0], 'y':graph[:,1]}

    return True