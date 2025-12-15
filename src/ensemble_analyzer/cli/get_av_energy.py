import argparse
import json
import os
import sys
import numpy as np
from pathlib import Path

from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer._managers.checkpoint_manager import CheckpointManager
from ensemble_analyzer._managers.protocol_manager import ProtocolManager
from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer.rrho import free_gibbs_energy
from ensemble_analyzer.constants import EH_TO_KCAL, R, CAL_TO_J
from ensemble_analyzer.title import title

def get_thermo_data(conf: Conformer, protocol_number: int, temp: float, mult: int, cut_off:float, alpha: int, pressure: float, linear:bool):
    """
    Retrieves Electronic Energy (E) and attempts to calculate ZPVE, H, G.
    Returns a tuple (E, E_ZPVE, H, G).
    If frequencies are missing, E_ZPVE, H, and G will be np.nan.
    """
    
    # 1. Retrieve Electronic Energy (E)
    if not conf.energies.__contains__(protocol_number):
        return np.nan, np.nan, np.nan, np.nan
    
    record_curr = conf.energies.__getitem__(protocol_number)
    E = record_curr.E
    
    # 2. Search for Frequencies (Current step or recursive fallback)
    freq = None
    for step in range(int(protocol_number), -1, -1):
        if conf.energies.__contains__(step):
            r = conf.energies.__getitem__(step)
            if r.Freq is not None and len(r.Freq) > 0:
                freq = r.Freq
                break
    
    # If no frequencies found, return only E
    if freq is None:
        return E, np.nan, np.nan, np.nan

    # 3. Necessary data for RRHO
    mw = conf.weight_mass
    B_vec = record_curr.B_vec
    
    # Fallback search for B_vec (Rotational Constants)
    if B_vec is None:
         for step in range(int(protocol_number), -1, -1):
            if conf.energies.__contains__(step):
                r = conf.energies.__getitem__(step)
                if r.B_vec is not None:
                    B_vec = r.B_vec
                    break
    
    # Dummy fallback to avoid crashes if B is missing but Freq exists (rare edge case)
    if B_vec is None:
        B_vec = np.array([1.0, 1.0, 1.0]) 

    # 4. Thermodynamic Calculation
    try:
        # G = Gibbs, h = Enthalpy (SCF + ZPVE + Therm_Corr)
        G, zpve, h, S = free_gibbs_energy(
            SCF=E, T=temp, freq=freq, mw=mw, B=B_vec, m=mult, cut_off=cut_off, alpha=alpha, P=pressure, linear=linear
        )
        return E, E + zpve, E + h, G
    except Exception:
        # In case of math errors 
        return E, np.nan, np.nan, np.nan

def calculate_population_vector(energies: np.ndarray, temp: float):
    """
    Calculates Boltzmann population (%) vector given an energy vector (Eh).
    Handles NaNs by returning NaNs in corresponding positions.
    """
    mask = ~np.isnan(energies)
    if not np.any(mask):
        return np.full(energies.shape, np.nan)
    
    valid_energies = energies[mask]
    
    # Delta E in kcal/mol relative to the valid minimum
    rel_energies = (valid_energies - valid_energies.min()) * EH_TO_KCAL
    
    # Boltzmann: exp(-dE / RT)
    exponent = -(rel_energies * 1000 * CAL_TO_J) / (R * temp)
    weights = np.exp(exponent)
    
    valid_pops = weights / weights.sum() * 100
    
    # Reconstruct full vector
    pops = np.full(energies.shape, np.nan)
    pops[mask] = valid_pops
    return pops

def calculate_weighted_average(energies: np.ndarray, pops: np.ndarray):
    """Calculates weighted average ignoring NaNs."""
    mask = ~np.isnan(energies) & ~np.isnan(pops)
    if not np.any(mask):
        return np.nan
    
    e_valid = energies[mask]
    p_valid = pops[mask]
    
    # p_valid is in %, divide by 100
    return np.sum(e_valid * (p_valid / 100.0))

def main():
    parser = argparse.ArgumentParser(description="Multi-Level Average Energy Analysis")
    parser.add_argument("-d", "--dir", default=".", help="Working directory")
    parser.add_argument("-T", "--temp", type=float, help="Temperature (K) for recalculation.")
    parser.add_argument("-o", "--output", default="average_energy_report.log", help="Output log file")

    # Thermo Parameters
    parser.add_argument("--cut-off", type=float, default=100.0, 
                        help="qRRHO cut-off frequency [cm-1]. Default: 100.0")
    parser.add_argument("--alpha", type=int, default=4, 
                        help="qRRHO damping factor alpha. Default: 4")
    parser.add_argument("--pressure", type=float, default=101.325, 
                        help="Pressure [kPa]. Default: 101.325")
    parser.add_argument('--linear', help='Define if molecules are linear', 
                        action='store_true')
    
    # Arguments for mathematical operations
    parser.add_argument("--sub", nargs=2, action='append', metavar=('P1', 'P2'), 
                        help="Subtraction: Avg(P1) - Avg(P2). Example: --sub 4 2")
    parser.add_argument("--add", nargs=2, action='append', metavar=('P1', 'P2'), 
                        help="Addition: Avg(P1) + Avg(P2). Example: --add 1 3")
    
    args = parser.parse_args()
    work_dir = Path(args.dir)
    
    # Setup Logger
    logger = create_logger(Path(args.output), debug=False)
    
    # Load Settings
    settings_path = work_dir / "settings.json"
    if not settings_path.exists():
        logger.critical(f"Settings file not found at {settings_path}")
        sys.exit(1)
        
    with open(settings_path, 'r') as f:
        settings = json.load(f)
    
    original_temp = settings.get("temperature", 298.15)
    target_temp = args.temp if args.temp is not None else original_temp
    

    logger.info(title)
    logger.info(f"Analysis Temperature: {target_temp} K")
    if abs(target_temp - original_temp) > 1e-3:
        logger.info("Performing thermodynamic recalculation due to temperature change.")

    # Load Data (Context Manager for directory switching)
    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        ckpt_mgr = CheckpointManager()
        conformers = ckpt_mgr.load()
        protocol_mgr = ProtocolManager( )
        protocols = protocol_mgr.load()
    finally:
        os.chdir(cwd)

    # Sort protocols numerically
    protocols.sort(key=lambda x: int(x.number))
    
    final_summary_rows = []
    
    # Dictionary to store raw averages (Hartree) for math operations
    # Format: { p_num: {'E': float, 'EZPVE': float, 'H': float, 'G': float} }
    protocol_averages = {}

    # --- MAIN PROTOCOL LOOP ---
    for proto in protocols:
        p_num = int(proto.number)
        
        # 1. Collect Conformer Data
        data_rows = []
        
        for c in conformers:
            # Check if calculation exists for this step and is not being deactivated
            if not c.energies.__contains__(p_num):
                continue
            if np.isnan(c.energies.__getitem__(p_num).Pop): 
                continue
                
            # Retrieve Energies (Recalculated on the fly)
            e_val, ezpve_val, h_val, g_val = get_thermo_data(c, p_num, target_temp, int(proto.mult), cut_off=args.cut_off, alpha=args.alpha, pressure=args.pressure, linear=args.linear)
            
            # If even E is missing, skip
            if np.isnan(e_val):
                continue
                
            # Global activity check
            if not c.active:
                continue

            data_rows.append({
                "conf_obj": c,
                "E": e_val,
                "E_ZPVE": ezpve_val,
                "H": h_val,
                "G": g_val
            })

        if not data_rows:
            logger.warning(f"Protocol {p_num}: No active conformers found.")
            continue

        # 2. Calculate Populations
        vec_E = np.array([d["E"] for d in data_rows])
        vec_EZPVE = np.array([d["E_ZPVE"] for d in data_rows])
        vec_H = np.array([d["H"] for d in data_rows])
        vec_G = np.array([d["G"] for d in data_rows])
        
        pop_E = calculate_population_vector(vec_E, target_temp)
        pop_EZPVE = calculate_population_vector(vec_EZPVE, target_temp)
        pop_H = calculate_population_vector(vec_H, target_temp)
        pop_G = calculate_population_vector(vec_G, target_temp)

        # 3. Build Detailed Table
        table_rows = []
        for i, d in enumerate(data_rows):
            # Formatter for Energies (10 decimal places)
            fmt = lambda x: f"{x:.10f}" if not np.isnan(x) else "  ---  "
            # Formatter for Population (2 decimal places)
            fmt_pop = lambda x: f"{x:5.2f}" if not np.isnan(x) else " --- "
            
            comment = getattr(proto, 'comment', '')

            row = [
                f"{d['conf_obj'].number}",
                fmt(vec_E[i]),      fmt_pop(pop_E[i]),
                fmt(vec_EZPVE[i]),  fmt_pop(pop_EZPVE[i]),
                fmt(vec_H[i]),      fmt_pop(pop_H[i]),
                fmt(vec_G[i]),      fmt_pop(pop_G[i]),
            ]
            table_rows.append(row)
        
        headers = [
            "Conf", 
            "E [Eh]", "Pop(E)%", 
            "E+ZPVE", "Pop(EZ)%", 
            "H [Eh]", "Pop(H)%", 
            "G [Eh]", "Pop(G)%"
        ]
        
        logger.table(
            title=f"Protocol {p_num} Analysis @ {target_temp}K",
            headers=headers,
            data=table_rows,
            char="-"
        )

        # 4. Calculate Weighted Averages
        av_E = calculate_weighted_average(vec_E, pop_E)
        av_EZPVE = calculate_weighted_average(vec_EZPVE, pop_EZPVE)
        av_H = calculate_weighted_average(vec_H, pop_H)
        av_G = calculate_weighted_average(vec_G, pop_G)
        
        # Store for Operations
        protocol_averages[p_num] = {
            "E": av_E, "EZPVE": av_EZPVE, "H": av_H, "G": av_G,
        }

        # 5. Prepare Final Summary Row
        fmt_av = lambda x: f"{x:.10f}" if not np.isnan(x) else "---"
        
        final_summary_rows.append([
            f"{p_num}",
            f"{proto.functional}/{proto.basis}",
            comment,
            fmt_av(av_E),
            fmt_av(av_EZPVE),
            fmt_av(av_H),
            fmt_av(av_G),
            len(vec_E)
        ])

    # --- FINAL SUMMARY TABLE ---
    summary_headers = [
        "Prot.", "Level", "Comment", 
        "E_av [Eh]", "(E+ZPVE)_av", "H_av [Eh]", "G_av [Eh]", "N Conf."
    ]
    
    logger.table(
        title=f"Ensemble Average Energies Summary (Hartree) @ {target_temp} K",
        headers=summary_headers,
        data=final_summary_rows,
        char="="
    )

    # --- MATH OPERATIONS SECTION ---
    if args.sub or args.add:
        ops_rows = []
        
        def perform_op(p1_str, p2_str, op_type):
            try:
                p1, p2 = int(p1_str), int(p2_str)
            except ValueError:
                return [f"{op_type} {p1_str} {p2_str}", "Error: Invalid ID", "", "", ""]

            if p1 not in protocol_averages or p2 not in protocol_averages:
                return [f"{op_type} {p1} {p2}", "Error: Missing Data", "", "", ""]

            v1 = protocol_averages[p1]
            v2 = protocol_averages[p2]
            
            # Conversion factor: Output differences in kcal/mol for readability
            factor = EH_TO_KCAL 
            
            row = [f"Prot {p1} {op_type} {p2}"]
            for key in ["E", "EZPVE", "H", "G"]:
                if np.isnan(v1[key]) or np.isnan(v2[key]):
                    row.append("NaN")
                else:
                    val = (v1[key] - v2[key]) if op_type == "-" else (v1[key] + v2[key])
                    # Difference formatted to 2 decimal places (kcal/mol)
                    row.append(f"{val * factor:.2f}")
            return row

        # Perform Subtractions
        if args.sub:
            for p1, p2 in args.sub:
                ops_rows.append(perform_op(p1, p2, "-"))
        
        # Perform Additions
        if args.add:
            for p1, p2 in args.add:
                ops_rows.append(perform_op(p1, p2, "+"))

        if ops_rows:
            logger.table(
                title=f"Calculated Differences/Sums [kcal/mol] @ {target_temp} K",
                headers=["Operation", "∆E", "∆(E+ZPVE)", "∆H", "∆G"],
                data=ops_rows,
                char="*"
            )

if __name__ == "__main__":
    main()