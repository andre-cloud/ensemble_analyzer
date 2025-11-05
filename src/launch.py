from src.conformer import Conformer
from src.ioFile import read_ensemble, save_snapshot
from src.logger import create_log, ordinal, DEBUG
from src.parser_arguments import parser_arguments
from src.parser_parameter import get_conf_parameters
from src.IOsystem import SerialiseEncoder, move_files
from src.protocol import Protocol, load_protocol
from src.pruning import calculate_rel_energies, check_ensemble
from src.graph import main_graph, Compared
from src.regex_parsing import regex_parsing

from src.clustering import perform_PCA, get_ensemble
from src.title import title

from src.constants import R, MAX_TRY, EH_TO_KCAL

import time
import json
import datetime
import logging
from tabulate import tabulate
import os
from typing import Union, List
import numpy as np


def launch(idx, conf, protocol, cpu, log, temp, ensemble, try_num: int = 1) -> None:
    """
    Run the calculation for each conformer

    :param idx: index of the calculation
    :type idx: int
    :param conf: conformer instance
    :type conf: Conformer
    :param protocol: protocol instance
    :type protocol: Protocol
    :param cpu: number of cpu to allocate
    :type cpu: int
    :param log: instance
    :type log: logger
    :param temp: temperature [K]
    :type temp: float
    :param ensemble: whole ensemble list
    :type ensemble: list

    :return: None
    """

    log.info(
        f"{idx}. Running {ordinal(int(protocol.number))} PROTOCOL -> CONF{conf.number}"
    )

    # Run actual calculation
    calc, label = protocol.get_calculator(cpu=cpu, conf=conf)
    atoms = conf.get_ase_atoms(calc)
    st = time.perf_counter()
    try:
        atoms.get_potential_energy()
    except Exception as e:
        pass
    end = time.perf_counter()
    ####

    move_files(conf, protocol, label)

    output_file = os.path.join(
        os.getcwd(),
        conf.folder,
        f"protocol_{protocol.number}",
        f'{conf.number}_p{protocol.number}_{label}.{regex_parsing[protocol.calculator]["ext"]}',
    )

    if not get_conf_parameters(
        conf=conf,
        number=protocol.number,
        output=output_file,
        p=protocol,
        time=end - st,
        temp=temp,
        log=log,
    ):

        # È da capire PERCHÉ il conformero non ha i parametri...
        # Qui i parametri sono già stati divisi e parsati, quindi il calcolo ha finito. Questo è il punto in si può fare il displacement sulla frequenza negativa (da capire la se TS oppure no...) o altre manipolazioni se il conto non ha trovato il minimo/sella.

        if try_num <= MAX_TRY:

            log.error(
                f"ERROR: During calculation of CONF_{conf.number} a server error occur and the energy could not be parsed; re-running protocol {protocol.number} on the same conformer for the {ordinal(try_num)} time"
            )
            time.sleep(10)
            launch(idx, conf, protocol, cpu, log, temp, ensemble, try_num=try_num + 1)
        else:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nMax number of re-run ({MAX_TRY}) executed for CONF_{conf.number}.\n{'='*20}\nExiting\n{'='*20}"
            )
            raise RuntimeError(
                f"Max number of re-run ({MAX_TRY}) executed for CONF_{conf.number}. Exiting"
            )

    json.dump(
        {i.number: i.__dict__ for i in ensemble},
        open("checkpoint.json", "w"),
        indent=4,
        cls=SerialiseEncoder,
    )


# exclude_enantiomers
def run_protocol(
    conformers: List[Conformer],
    p: Protocol,
    temperature: float,
    cpu: int,
    log,
    include_H: bool,
) -> None:
    """
    Run the protocol for each conformer

    :param conformers: whole ensemble list
    :type conformers: list
    :param p: protocol information
    :type p: Protocol
    :param temperature: temperature [K]
    :type temperature: float
    :param cpu: cpu to allocate
    :type cpu: int
    :param log: logger instance
    :type log: logger

    :return: None
    """

    log.info(f"STARTING PROTOCOL {p.number}")
    log.info(
        f"\nActive conformers for this phase: {len([i for i in conformers if i.active])}\n"
    )

    # No enumerate cause conformers are both active and inactive

    count = 1
    for i in conformers:
        if not i.active:
            continue
        if i.energies.get(str(p.number)):
            continue
        launch(count, i, p, cpu, log, temperature, conformers)
        count += 1

    conformers = sorted(conformers)

    log.info("\nEnded Calculations\n")
    log.info(
        f"\nTotal elapsed time for protocol {ordinal(int(p.number))}: "
        + str(
            datetime.timedelta(
                seconds=sum([i._last_energy["time"] for i in conformers if i.active])
            )
        )
    )

    if DEBUG:
        if p.cluster:
            perform_PCA(
                confs=[i for i in conformers if i.active],
                ncluster=p.cluster if type(p.cluster) is int else None,
                fname=f"PCA_before_pruning_protocol_{p.number}.png",
                title=f"PCA before pruning protocol {p.number}",
                log=log,
                set=False,
                include_H=include_H,
            )
        create_summary("Summary", conformers, p, log)

    log.debug("Start Pruning")

    # , exclude_enantiomers)
    conformers = check_ensemble(conformers, p, log, include_H)
    conformers = sort_conformers_by_energy(conformers, temperature)

    save_snapshot(f"ensemble_after_{p.number}.xyz", conformers, log)

    if p.cluster:
        perform_PCA(
            confs=[i for i in conformers if i.active],
            ncluster=p.cluster if type(p.cluster) is int else None,
            fname=f"PCA_after_pruning_protocol_{p.number}.png",
            title=f"PCA after pruning protocol {p.number}",
            log=log,
            include_H=include_H,
            set_=True,
        )
    if type(p.cluster) is int:
        conformers = get_ensemble(conformers)

    create_summary("Summary After Pruning", conformers, p, log)

    calc_average_ensemble(conformers, p.number, temperature, log)

    log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    return None


def calc_average_ensemble(conformers: list, number, T, log) -> None:

    # ∆E, ∆(E+ZPVE), ∆H, ∆G
    CONFS = [i for i in conformers if i.active]

    dE = np.array([i.energies[str(number)]["E"] for i in CONFS])
    dE_ZPVE = np.array(
        [
            i.energies[str(number)]["E"] + i.energies[str(number)]["zpve"] * EH_TO_KCAL
            for i in CONFS
        ]
    )
    dH = np.array(
        [
            i.energies[str(number)]["E"] + i.energies[str(number)]["H"] * EH_TO_KCAL
            for i in CONFS
        ]
    )
    dG = np.array([i.energies[str(number)]["G"] for i in CONFS])

    # Boltzmann populations
    dE_boltz = boltz(dE, T)
    dE_ZPVE_boltz = boltz(dE_ZPVE, T)
    dH_boltz = boltz(dH, T)
    dG_boltz = boltz(dG, T)

    averages = [
        float(np.sum(dE * dE_boltz)) / EH_TO_KCAL,
        float(np.sum(dE_ZPVE * dE_ZPVE_boltz)) / EH_TO_KCAL,
        float(np.sum(dH * dH_boltz)) / EH_TO_KCAL,
        float(np.sum(dG * dG_boltz)) / EH_TO_KCAL,
    ]

    rows = [
        [
            f"Conf {i}",
            dE[idx] / EH_TO_KCAL,
            f"{round(dE_boltz[idx]*100, 2)}",
            dE_ZPVE[idx] / EH_TO_KCAL,
            f"{round(dE_ZPVE_boltz[idx]*100, 2)}",
            dH[idx] / EH_TO_KCAL,
            f"{round(dH_boltz[idx]*100, 2)}",
            dG[idx] / EH_TO_KCAL,
            f"{round(dG_boltz[idx]*100, 2)}",
        ]
        for idx, i in enumerate(CONFS)
    ]
    log.info("Energetic Summary of the active conformers")
    log.info(
        tabulate(
            rows,
            headers=[
                "Conformer",
                "∆E [Eh]",
                "Boltzamnn Pop. on ∆E",
                "∆(E+ZPVE) [Eh]",
                "Boltzamnn Pop. on ∆(E+ZPVE)",
                "∆H [Eh]",
                "Boltzamnn Pop. on ∆H",
                "∆G [Eh]",
                "Boltzamnn Pop. on ∆G",
            ],
            floatfmt=".10f",
        )
    )
    log.info("\n")
    log.info("Ensemble Avarages")
    log.info(
        tabulate(
            (averages,),
            headers=["E_av [Eh]", "E+ZPVE_av [Eh]", "H_av [Eh]", "G_av [Eh]"],
            floatfmt=".10f",
        )
    )
    log.info("\n\n")

    return


def boltz(energy: np.ndarray, T):
    ens = np.array(energy)
    ens -= min(ens)
    bolz = np.exp((-ens * 4186) / (R * T))
    return bolz / np.sum(bolz)


def last_protocol_completed(conf, idx: int) -> bool:
    """
    Getting if all the conformers have been calculated until the idx-th protocol

    :param conf: whole ensemble list
    :type conf: conf
    :param idx: index of the last protocol executed
    :type idx: idx

    :return: return if the protocol selected is completed
    :rtype: bool
    """

    tmp = []
    for i in conf:
        if i.energies.get(int(idx)) is not None and i.active:
            tmp.append(i)

    return len([tmp]) == 0


def create_summary(title, conformers, protocol, log):
    """
    Create the summary of a ensemble information

    :param title: title of the summary
    :type title: str
    :param conformers: whole ensemble list
    :type conformers: list
    :param log: logger instance
    :type log: logger

    :return: None
    """

    headers = [
        "Conformers",
        "E [Eh]",
        "G-E [Eh]",
        "G [Eh]",
        "B [cm-1]",
        "∆G [kcal/mol]",
        "Pop [%]",
        "Elap. time [sec]",
        "# Cluster",
    ] + [i for i in list(protocol.verbal_internals())]


    log.info(title)
    log.info("")
    log.info(
        tabulate(
            [i.create_log(protocol.monitor_internals) for i in conformers if i.active],
            headers=headers,
            floatfmt=".5f",
        )
    )
    log.info("\n\n")

    return None


def start_calculation(
    conformers,
    protocol,
    cpu: int,
    temperature: float,
    start_from: int,
    log,
    final_lambda=800.0,
    definition=4,
    FWHM: Union[None, float, List[float]] = None,
    shift: Union[None, float, List[float]] = None,
    invert=False,
    include_H=True,
    # exclude_enantiomers=False,
) -> None:
    """
    Main calculation loop

    :param conformers: whole ensemble list
    :type conformers: list
    :param protocol: whole protocol steps
    :type protocol: list
    :param cpu: cpu allocated
    :type cpu: int
    :param temperature: temperature [K]
    :type temperature: float
    :param start_from: index of the last protocol executed
    :type start_from: int
    :param log: logger instance
    :type log: logger
    :param include_H: include the hydrogen atoms in the PCA and RSMD
    :type include_H: bool
    :param exclude_enantiomers: exclude enantiomeric conformers
    :type exclude_enantiomers: bool

    :return: None
    """
    if start_from != 0:
        if last_protocol_completed(conformers, start_from):
            conformers = check_ensemble(conformers, protocol[start_from], log)
            create_summary("Summary", conformers, log)

    for p in protocol[start_from:]:
        with open("last_protocol", "w") as f:
            f.write(str(p.number))

        run_protocol(
            conformers, p, temperature, cpu, log, include_H
        )  # , exclude_enantiomers)

        log.debug("Creating graph")

        # TODO: incorporate FWHM and shift from settings as bounds
        main_graph(conformers, p, log, invert=invert, read_pop=p.read_population)

    # sort the final ensemble
    c_ = sort_conformers_by_energy(conformers, temperature)
    save_snapshot("final_ensemble.xyz", c_, log)

    create_summary("Final Summary", c_, p, log)

    log.info("\n\n\n\n")
    log.info(
        f"Final ensemble has {len([i for i in conformers if i.active])} conformers"
    )

    t = 0
    for i in conformers:
        for j in i.energies:
            t += i.energies[j]["time"]
    log.info(f"Total elapsed time: {datetime.timedelta(seconds=t)}")

    log.info("Ensemble refined correctly!\n\n")
    log.info(f'{"="*15}\nCALCULATIONS ENDED\n{"="*15}\n\n')

    # save the final ensemble to a checkpoint file
    json.dump(
        {i.number: i.__dict__ for i in conformers},
        open("checkpoint.json", "w"),
        indent=4,
        cls=SerialiseEncoder,
    )

    return None


def sort_conformers_by_energy(conformers, temperature):
    calculate_rel_energies(conformers, temperature)
    c_ = sorted(conformers)
    return c_


def restart() -> tuple:
    """
    Reload ensemble, protocol and setting from previous calculation

    :return: the ensemble list, the protocol list and the last completed protocol
    :rtype: list, list, int
    """
    confs = json.load(open("checkpoint.json"))
    ensemble = [Conformer.load_raw(confs[i]) for i in confs]

    p = json.load(open("protocol_dump.json"))
    protocol = [Protocol(**p[i]) for i in p]

    with open("last_protocol") as f:
        start_from = int(f.readlines()[0])

    return ensemble, protocol, start_from


def create_protocol(p, log) -> list:
    """
    Create the steps for the protocol to be executed

    :param p: JSON read file of the protocol
    :type p: dict
    :param thrs: JSON read file of the thresholds
    :type thrs: dict
    :param log: logger instance
    :type log: logger

    :return: protocol steps
    :rtype: list
    """

    protocol = []

    log.debug("Loading Protocol\n")
    last_prot_with_freq = None
    for idx, d in p.items():
        func = d.get("functional", None)
        add_input = d.get("add_input", "")
        graph = d.get("graph", False)
        freq = d.get("freq", False)

        if not graph and (freq or "freq" in add_input.lower()):
            last_prot_with_freq = int(idx)

        check_protocol(log, func, graph, freq, add_input, idx, last_prot_with_freq)

        protocol.append(Protocol(number=idx, **d))

    log.info(
        ""
        "\n".join(
            (
                f"{i.number}: {str(i)} - {i.calculation_level}\n {i.thr}"
                for i in protocol
            )
        )
        + "\n"
    )
    return protocol


def check_protocol(
    log, func, graph, freq, add_input, idx, last_prot_with_freq=None
) -> None:
    """
    Sanity check of the input settings

    :param log: logger instance
    :type log: logger
    :param func: the functional name
    :type func: str
    :param graph: calculate the TD-DFT graphs
    :type graph: bool
    :param add_input: additional input to put in the input of the calculation
    :type add_input: str
    :param idx: index of the protocol
    :type idx: int
    :param last_prot_with_freq: last index with frequency calculation, default is None
    :type last_prot_with_freq: int, optional

    :return: None

    :raise IOError:
        There is an error in the input file with the definition of the functional. See the output file.
    """

    if not func:
        log.critical(
            f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNCTIONAL key must be passed in order to calculate energy.\nDFT functional, HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\nExiting\n{'='*20}\n"
        )
        raise IOError(
            "There is an error in the input file with the definition of the functional. See the output file."
        )

    if graph:
        if not add_input and not freq:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nADD_INPUT must be set so the proper calculation (TD-DFT or CASSCF/RASSCF) to simulate the electronic spectra (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\n Exiting\n{'='*20}\n"
            )
            raise IOError(
                "There is an error in the input file with the definition of the functional. See the output file."
            )
        if not isinstance(last_prot_with_freq, int):
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nElectrical spectra requires Boltzmann population over ∆G. In the specified protocol there is NO frequency calculation turned on (Problem at {ordinal(int(idx))} protocol definition).\n{'='*20}\nExiting\n{'='*20}"
            )
            raise IOError(
                "There is an error in the input file with the definition of the functional. See the output file."
            )

    return None


def main():
    """
    Main script loop
    """

    args = parser_arguments()

    # Trying to reload the damped setting from a previously calculation. Else damp the settings

    if os.path.exists(os.path.join(os.getcwd(), "settings.json")):
        settings = json.load(open("settings.json"))
    else:
        settings = {
            "output": args.output,
            "cpu": args.cpu,
            "temperature": args.temperature,
            "final_lambda": args.final_lambda,
            "definition": args.definition,
            "fwhm": args.fwhm,
            "shift": args.shift,
            "invert": args.invert,
            "include_H": not args.exclude_H,
            # "exclude_enantiomers" : args.exclude_enantiomers,
        }
        json.dump(settings, open("settings.json", "w"), indent=4)

    # create the setting dictionary
    output = (
        settings.get("output", args.output)
        if not args.restart
        else ".".join(settings.get("output", args.output).split(".")[:-1])
        + "_restart.out"
    )

    # initiate the log
    log = create_log(output)
    # deactivate the log of matplotlib
    logging.getLogger("matplotlib").disabled = False

    log.info(title)

    if args.restart:
        # reload the previous information from checkpoint file
        conformers, protocol, start_from = restart()

    else:
        protocol = create_protocol(load_protocol(args.protocol), log)
        start_from = protocol[0].number
        json.dump(
            {i.number: i.__dict__ for i in protocol},
            open("protocol_dump.json", "w"),
            indent=4,
            cls=SerialiseEncoder,
        )
        conformers = read_ensemble(args.ensemble, log)
        if len(conformers) > 10:
            perform_PCA(
                conformers, None, "initial_pca", "PCA analysis of Conf Search", log
            )

    # start the loop
    start_calculation(
        conformers=conformers,
        protocol=protocol,
        start_from=int(start_from),
        log=log,
        cpu=settings.get("cpu", args.cpu),
        temperature=settings.get("temperature", args.temperature),
        final_lambda=settings.get("final_lambda", 800),
        definition=settings.get("definition", 4),
        FWHM=settings.get("fwhm", None),
        shift=settings.get("shift", None),
        invert=settings.get("invert", False),
        include_H=settings.get("include_H", True),
        # exclude_enantiomers=settings.get("exclude_enantiomers", False),
    )

    for j in ["IR", "VCD", "UV", "ECD"]:
        g = Compared(protocol, graph_type=j)
        g.save_graph()


if __name__ == "__main__":
    conformers, protocol, start_from = restart()
    print(f"Final ensemble has {len(conformers)} conformers")
    t = 0
    for i in conformers:
        for j in i.energies:
            t += i.energies[j]["time"]
    print(f"Total elapsed time: {datetime.timedelta(seconds=t)}")
