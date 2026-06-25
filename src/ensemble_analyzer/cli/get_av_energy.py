import argparse
import json
import os
import sys
from pathlib import Path
from typing import Tuple


def get_thermo_data(conf, protocol_number, temp, mult, cut_off, alpha, pressure, linear) -> Tuple[float, ...]:
    import numpy as np
    from ensemble_analyzer.rrho import free_gibbs_energy

    if protocol_number not in conf.energies:
        return np.nan, np.nan, np.nan, np.nan

    record_curr = conf.energies[protocol_number]
    E = record_curr.E

    freq = conf.energies.get_last_freq(int(protocol_number))
    if len(freq) == 0:
        return E, np.nan, np.nan, np.nan

    mw = conf.weight_mass
    B_vec = conf.energies.get_last_bvec(int(protocol_number))
    if B_vec is None:
        B_vec = np.array([1.0, 1.0, 1.0])

    try:
        G, zpve, h, S = free_gibbs_energy(
            SCF=E, T=temp, freq=freq, mw=mw, B=B_vec, m=mult,
            cut_off=cut_off, alpha=alpha, P=pressure, linear=linear,
        )
        return E, E + zpve, E + h, G
    except Exception:
        return E, np.nan, np.nan, np.nan


def calculate_population_vector(energies, temp):
    import numpy as np
    from ensemble_analyzer.constants import boltzmann_distribution

    mask = ~np.isnan(energies)
    if not np.any(mask):
        return np.full(energies.shape, np.nan)
    _, valid_pops = boltzmann_distribution(energies[mask], temp)
    pops = np.full(energies.shape, np.nan)
    pops[mask] = valid_pops * 100
    return pops


def calculate_weighted_average(energies, pops):
    import numpy as np

    mask = ~np.isnan(energies) & ~np.isnan(pops)
    if not np.any(mask):
        return np.nan

    e_valid = energies[mask]
    p_valid = pops[mask]

    return np.sum(e_valid * (p_valid / 100.0))


def main() -> None:
    parser = argparse.ArgumentParser(description="Multi-Level Average Energy Analysis")
    parser.add_argument("-d", "--dir", default=".", help="Working directory")
    parser.add_argument("-T", "--temp", type=float, help="Temperature (K) for recalculation.")
    parser.add_argument("-o", "--output", default="average_energy_report.log", help="Output log file")

    parser.add_argument("--cut-off", type=float, default=100.0,
                        help="qRRHO cut-off frequency [cm-1]. Default: 100.0")
    parser.add_argument("--alpha", type=int, default=4,
                        help="qRRHO damping factor alpha. Default: 4")
    parser.add_argument("--pressure", type=float, default=101.325,
                        help="Pressure [kPa]. Default: 101.325")
    parser.add_argument('--linear', help='Define if molecules are linear',
                        action='store_true')

    parser.add_argument("--sub", nargs=2, action='append', metavar=('P1', 'P2'),
                        help="Subtraction: Avg(P1) - Avg(P2). Example: --sub 4 2")
    parser.add_argument("--add", nargs=2, action='append', metavar=('P1', 'P2'),
                        help="Addition: Avg(P1) + Avg(P2). Example: --add 1 3")

    parser.add_argument("--validate", nargs=4, action='append',
                        metavar=('Protocollo', 'Pattern', 'Value', 'Thr'),
                        help="Validate conformer output post-hoc: "
                             "--validate Protocollo 'regex' expected threshold")

    args = parser.parse_args()

    from collections import defaultdict
    from ensemble_analyzer._logger.create_log import create_logger
    from ensemble_analyzer.ensemble_io import load_workflow_data
    from ensemble_analyzer.protocol.protocol import sort_protocols
    from ensemble_analyzer._title import title
    from ensemble_analyzer.validators import validate_line
    from ensemble_analyzer.constants import regex_parsing

    validators_by_proto: dict[str, list] = defaultdict(list)
    if args.validate:
        for proto, pattern, val_str, thr_str in args.validate:
            validators_by_proto[proto].append((pattern, float(val_str), float(thr_str)))

    work_dir = Path(args.dir)

    logger = create_logger(Path(args.output), debug=False)

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

    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        conformers, protocols = load_workflow_data()
    finally:
        os.chdir(cwd)

    protocols = sort_protocols(protocols)

    final_summary_rows = []
    protocol_averages = {}

    for proto in protocols:
        p_num = int(proto.number)

        data_rows = []
        failed_validations = []

        for c in conformers:
            import numpy as np
            if p_num not in c.energies:
                continue
            record = c.energies[p_num]
            if np.isnan(record.Pop):
                continue

            proto_validators = validators_by_proto.get(str(p_num))
            if proto_validators:
                calc_name = record.calculator.lower()
                ext = regex_parsing.get(calc_name, {}).get("ext")
                if ext is None:
                    logger.warning(f"ML calculator {calc_name}: post-hoc validate skipped for conf {c.number}")
                else:
                    proto_dir = Path(c.folder) / f"protocol_{p_num}"
                    matches = list(proto_dir.glob(f"{c.number}_p{p_num}_*.{ext}"))
                    if not matches:
                        logger.warning(f"Output not found for conf {c.number}, proto {p_num}: {proto_dir}")
                        continue
                    text = matches[0].read_text()
                    ok = all(
                        validate_line(text, pattern, expected, threshold)
                        for pattern, expected, threshold in proto_validators
                    )
                    if not ok:
                        failed_validations.append(c.number)
                        continue

            e_val, ezpve_val, h_val, g_val = get_thermo_data(
                c, p_num, target_temp, int(proto.mult),
                cut_off=args.cut_off, alpha=args.alpha,
                pressure=args.pressure, linear=args.linear,
            )

            if np.isnan(e_val):
                continue

            data_rows.append({
                "conf_obj": c,
                "E": e_val,
                "E_ZPVE": ezpve_val,
                "H": h_val,
                "G": g_val,
            })

        if not data_rows:
            logger.warning(f"Protocol {p_num}: No active conformers found.")
            continue

        import numpy as np
        vec_E = np.array([d["E"] for d in data_rows])
        vec_EZPVE = np.array([d["E_ZPVE"] for d in data_rows])
        vec_H = np.array([d["H"] for d in data_rows])
        vec_G = np.array([d["G"] for d in data_rows])

        pop_E = calculate_population_vector(vec_E, target_temp)
        pop_EZPVE = calculate_population_vector(vec_EZPVE, target_temp)
        pop_H = calculate_population_vector(vec_H, target_temp)
        pop_G = calculate_population_vector(vec_G, target_temp)

        table_rows = []
        for i, d in enumerate(data_rows):
            fmt = lambda x: f"{x:.10f}" if not np.isnan(x) else "  ---  "
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
            "G [Eh]", "Pop(G)%",
        ]

        logger.table(
            title=f"Protocol {p_num} Analysis @ {target_temp}K",
            headers=headers,
            data=table_rows,
            char="-",
        )

        if failed_validations:
            for conf_id in failed_validations:
                logger.warning(f"Validation failed for conf {conf_id}, proto {p_num}")

        av_E = calculate_weighted_average(vec_E, pop_E)
        av_EZPVE = calculate_weighted_average(vec_EZPVE, pop_EZPVE)
        av_H = calculate_weighted_average(vec_H, pop_H)
        av_G = calculate_weighted_average(vec_G, pop_G)

        protocol_averages[p_num] = {
            "E": av_E, "EZPVE": av_EZPVE, "H": av_H, "G": av_G,
        }

        fmt_av = lambda x: f"{x:.10f}" if not np.isnan(x) else "---"

        final_summary_rows.append([
            f"{p_num}",
            f"{proto.functional}/{proto.basis}",
            comment,
            fmt_av(av_E),
            fmt_av(av_EZPVE),
            fmt_av(av_H),
            fmt_av(av_G),
            len(vec_E),
        ])

    summary_headers = [
        "Prot.", "Level", "Comment",
        "E_av [Eh]", "(E+ZPVE)_av", "H_av [Eh]", "G_av [Eh]", "N Conf.",
    ]

    logger.table(
        title=f"Ensemble Average Energies Summary (Hartree) @ {target_temp} K",
        headers=summary_headers,
        data=final_summary_rows,
        char="=",
    )

    if args.sub or args.add:
        ops_rows = []

        def perform_op(p1_str, p2_str, op_type):
            import numpy as np
            from ensemble_analyzer.constants import EH_TO_KCAL

            try:
                p1, p2 = int(p1_str), int(p2_str)
            except ValueError:
                return [f"{op_type} {p1_str} {p2_str}", "Error: Invalid ID", "", "", ""]

            if p1 not in protocol_averages or p2 not in protocol_averages:
                return [f"{op_type} {p1} {p2}", "Error: Missing Data", "", "", ""]

            v1 = protocol_averages[p1]
            v2 = protocol_averages[p2]

            factor = EH_TO_KCAL

            row = [f"Prot {p1} {op_type} {p2}"]
            for key in ["E", "EZPVE", "H", "G"]:
                if np.isnan(v1[key]) or np.isnan(v2[key]):
                    row.append("NaN")
                else:
                    val = (v1[key] - v2[key]) if op_type == "-" else (v1[key] + v2[key])
                    row.append(f"{val * factor:.2f}")
            return row

        if args.sub:
            for p1, p2 in args.sub:
                ops_rows.append(perform_op(p1, p2, "-"))

        if args.add:
            for p1, p2 in args.add:
                ops_rows.append(perform_op(p1, p2, "+"))

        if ops_rows:
            logger.table(
                title=f"Calculated Differences/Sums [kcal/mol] @ {target_temp} K",
                headers=["Operation", "∆E", "∆(E+ZPVE)", "∆H", "∆G"],
                data=ops_rows,
                char="*",
            )


if __name__ == "__main__":
    main()
