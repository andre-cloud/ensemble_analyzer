try:
    from src.conformer import Conformer
    from src.ioFile import read_ensemble, save_snapshot
    from src.logger import create_log, ordinal
    from src.parser_arguments import parser_arguments
    from src.parser_parameter import get_conf_parameters
    from src.IOsystem import SerialiseEncoder
    from src.protocol import Protocol, load_protocol
    from src.pruning import calculate_rel_energies, check_ensemble
    from src.grapher import Graph, plot_conv_graph
    from src.clustering import perform_PCA
    from src.title import title
    from src.calculations import optimize, calc_freq, single_point
except ImportError as e:  # pragma: no cover
    print(e)
    from conformer import Conformer
    from ioFile import read_ensemble, save_snapshot
    from logger import create_log, ordinal
    from parser_arguments import parser_arguments
    from parser_parameter import get_conf_parameters
    from IOsystem import SerialiseEncoder
    from protocol import Protocol, load_protocol
    from pruning import calculate_rel_energies, check_ensemble
    from grapher import Graph
    from clustering import perform_PCA
    from title import title
    from calculations import optimize, calc_freq, single_point

import ase
import time
import json
import datetime
import logging
from tabulate import tabulate
import os
from typing import Union


MAX_TRY = 5


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
    st = time.perf_counter()

    if protocol.opt:
        atoms, label = optimize(conf, protocol, cpu=cpu, log=log)

    if protocol.freq:
        atoms, label = calc_freq(conf, protocol, cpu=cpu, log=log)

    if not (protocol.opt or protocol.freq):
        atoms, label = single_point(conf, protocol, cpu=cpu, log=log)

    end = time.perf_counter()

    os.rename(f"{label}.out", f"{conf.folder}/protocol_{protocol.number}.out")
    if protocol.freq:
        os.rename(f"{label}.hess", f"{conf.folder}/protocol_{protocol.number}.hess")
    os.remove(f"{label}.gbw")

    if not get_conf_parameters(conf, protocol.number, protocol, end - st, temp, log):
        if try_num <= MAX_TRY:
            log.error(
                f"ERROR: During calculation of CONF_{conf.number} a server error occur and the energy could not be parsed; re-running protocol {protocol.number} on the same conformer for the {ordinal(try_num)} time"
            )
            time.sleep(10)
            launch(idx, conf, protocol, cpu, log, temp, ensemble, try_num=try_num + 1)
        else:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nMax number of re-run ({MAX_TRY}) executed for CONF_{conf.number}.{'='*20}\nExiting\n{'='*20}"
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


def run_protocol(conformers, p, temperature, cpu, log) -> None:
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
        "\nTotal elapsed time: "
        + str(
            datetime.timedelta(
                seconds=sum([i._last_energy["time"] for i in conformers if i.active])
            )
        )
    )

    perform_PCA(
        [i for i in conformers if i.active],
        5,
        f"PCA_before_pruning_protocol_{p.number}.png",
        f"PCA before pruning protocol {p.number}",
        log,
    )

    create_summary("Summary", conformers, log)
    conformers = check_ensemble(conformers, p, log)
    log.debug("Start Pruning")
    calculate_rel_energies(conformers, temperature)
    save_snapshot(f"ensemble_after_{p.number}.xyz", conformers, log)

    perform_PCA(
        [i for i in conformers if i.active],
        5,
        f"PCA_after_pruning_protocol_{p.number}.png",
        f"PCA after pruning protocol {p.number}",
        log,
    )

    calculate_rel_energies(conformers, temperature)
    create_summary("Summary After Pruning", conformers, log)

    log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')

    return None


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


def create_summary(title, conformers, log):
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

    log.info(title)
    log.info("")
    log.info(
        tabulate(
            [i.create_log() for i in conformers if i.active],
            headers=[
                "Conformers",
                "E[Eh]",
                "G[Eh]",
                "B[cm-1]",
                "E. Rel [kcal/mol]",
                "Pop [%]",
                "Elap. time [sec]",
            ],
            floatfmt=".6f",
        )
    )
    log.info("")

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
    FWHM: Union[None, float] = None,
    shift: Union[None, float] = None,
    invert=False,
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

    :return: None
    """
    if start_from != 0:
        if last_protocol_completed(conformers, start_from):
            conformers = check_ensemble(conformers, protocol[start_from], log)
            create_summary("Summary", conformers, log)

    for p in protocol[start_from:]:
        with open("last_protocol", "w") as f:
            f.write(str(p.number))
        run_protocol(conformers, p, temperature, cpu, log)
        if p.graph:
            Graph(
                confs=conformers,
                protocol=p,
                log=log,
                T=temperature,
                final_lambda=final_lambda,
                definition=definition,
                FWHM=FWHM,
                shift=shift,
                invert=invert,
            )

    # sort the final ensemble
    c_ = sorted(conformers)
    calculate_rel_energies(c_, temperature)
    save_snapshot("final_ensemble.xyz", c_, log)
    log.info(f'{"="*15}\nCALCULATIONS ENDED\n{"="*15}\n\n')
    create_summary("Final Summary", c_, log)

    return None


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

        check_protocol(log, func, graph, add_input, idx, last_prot_with_freq)

        if not graph and d.get("freq"):
            last_prot_with_freq = int(idx)

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


def check_protocol(log, func, graph, add_input, idx, last_prot_with_freq=None) -> None:
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
            f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNCTIONAL key must be passed in order to calculate energy.\nDFT functional or HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\nExiting\n{'='*20}\n"
        )
        raise IOError(
            "There is an error in the input file with the definition of the functional. See the output file."
        )

    if graph:
        if not add_input:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nADD_INPUT must be set so the proper calculation (TD-DFT or CASSCF/RASSCF) to simulate the electronic spectra (Problem at {ordinal(int(idx))} protocol definition)\n{'='*20}\n Exiting\n{'='*20}\n"
            )
            raise IOError(
                "There is an error in the input file with the definition of the functional. See the output file."
            )
        if not last_prot_with_freq:
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nElectrical spectra requires Boltzmann population over ∆G. In the specified protocol there is NO frequency calculation turned on.\n{'='*20}\nExiting\n{'='*20}"
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
        }
        json.dump(settings, open("settings.json", "w"), indent=4)

    # create the setting dictionary
    output = (
        settings.get("output", args.output)
        if not args.restart
        else ".".join(settings.get("output", args.output).split(".")[:-1])
        + "_restart.out"
    )
    cpu = settings.get("cpu", args.cpu)
    temperature = settings.get("temperature", args.temperature)
    final_lambda = settings.get("final_lambda", 800)
    definition = settings.get("definition", 4)
    fwhm = settings.get("fwhm", None)
    shift = settings.get("shift", None)
    invert = settings.get("invert", False)

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
        conformers = read_ensemble(args.ensemble, args.charge, args.multiplicity, log)
        if len(conformers) > 3:
            perform_PCA(
                conformers, 5, "initial_pca", "PCA analysis of Conf Search", log
            )

    # start the loop
    start_calculation(
        conformers=conformers,
        protocol=protocol,
        cpu=cpu,
        temperature=temperature,
        start_from=int(start_from),
        log=log,
        final_lambda=final_lambda,
        definition=definition,
        FWHM=fwhm,
        shift=shift,
        invert=invert,
    )

    # Plots eventual graphs
    idxs_graphs = [i.number for i in protocol if i.graph]
    if len(idxs_graphs) > 0:
        plot_conv_graph(idxs_graphs, protocol)

    log.info("Ensemble refined correcly!")
