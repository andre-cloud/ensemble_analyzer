import ase
import os
import re
from ase.constraints import FixAtoms
from ase.optimize import LBFGS
from ase.vibrations import Infrared

try:
    from src.IOsystem import tail
    from src.regex_parsing import regex_parsing
    from src.parser_parameter import get_opt_geometry
except ImportError:
    from IOsystem import tail
    from regex_parsing import regex_parsing
    from parser_parameter import get_opt_geometry
    
MAX_TRY = 5


def check_output(label, calc):
    """
    Check the output file

    :param label: label of the output
    :type label: str
    :param calc: name of the calculator
    :type calc: str
    :return: True if calculation is correctly finished
    :rtype: bool
    """

    file = f"{label}.{regex_parsing[calc]['ext']}"
    if os.path.exists(file):
        return regex_parsing[calc]["finish"] in tail(file, 5)

    return False


def optimize(conf, protocol, cpu: int, log, attempts=0):
    """
    Run an optimization calculation

    :param conf: the conformer instance to optimize
    :type conf: Conformer
    :param protocol: protocol instance to run
    :type protocol: Protocol
    :param cpu: number of cpu to run the optimization on
    :type cpu: int
    :param log: logger instance
    :type log: logging
    :param attempts: number of try of the calculation
    :type attempts: int
    :return: ase.atoms, label
    :rtype: tuple
    """

    calc, label = protocol.get_calculator(
        cpu=cpu, charge=conf.charge, mult=conf.mult, mode="opt"
    )

    atoms = conf.get_ase_atoms(calc)

    if protocol.constrains:
        c = FixAtoms(indices=list(protocol.constrains))
        atoms.set_constraint(c)

    opt = LBFGS(
        atoms,
        maxstep=protocol.maxstep,
        logfile=f"opt_p{protocol.number}_{conf.number}.log",
        trajectory=f"p{protocol.number}_{label}.trj",
    )
    opt.run(protocol.fmax, steps=1000)

    if not check_output(label, protocol.calculator):
        if attempts < MAX_TRY:
            optimize(conf, protocol, cpu, log, attempts + 1)
        else:
            with open(f"{label}.{regex_parsing[calc]['ext']}") as f:
                fl = f.read()
            log.error("\n".join(fl.splitlines()[-6:-3]))
            log.critical(
                f"\n{'='*20}\nCRITICAL ERROR\n{'='*20}\nSome sort of error have been encountered during the calculation of the calculator.\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise RuntimeError(
                "Some sort of error have been encountered during the calculation of the calculator."
            )

    set_last_geometry(conf, atoms.get_positions())

    os.rename(
        f"opt_p{protocol.number}_{conf.number}.log",
        f"{conf.folder}/opt_p{protocol.number}_{conf.number}.log",
    )
    os.rename(
        f"p{protocol.number}_{label}.trj",
        f"{conf.folder}/p{protocol.number}_{label}.trj",
    )

    return atoms, label


def calc_freq(conf, protocol, cpu: int, log, attempts=0):
    """
    Run an hessian calculation

    :param conf: the conformer instance to optimize
    :type conf: Conformer
    :param protocol: protocol instance to run
    :type protocol: Protocol
    :param cpu: number of cpu to run the optimization on
    :type cpu: int
    :param log: logger instance
    :type log: logging
    :param attempts: number of try of the calculation
    :type attempts: int
    :return: ase.atoms, label
    :rtype: tuple
    """

    calc, label = protocol.get_calculator(
        cpu=cpu, charge=conf.charge, mult=conf.mult, mode="freq"
    )

    atoms = conf.get_ase_atoms(calc)

    vib = Infrared(atoms)
    try:
        vib.run()
    except ase.calculators.calculator.PropertyNotImplementedError:
        pass
    except ase.calculators.calculator.ReadError:
        pass

    return atoms, label


def single_point(conf, protocol, cpu: int, log, attempts=0):
    """
    Run a single point calculation

    :param conf: the conformer instance to optimize
    :type conf: Conformer
    :param protocol: protocol instance to run
    :type protocol: Protocol
    :param cpu: number of cpu to run the optimization on
    :type cpu: int
    :param log: logger instance
    :type log: logging
    :param attempts: number of try of the calculation
    :type attempts: int
    :return: ase.atoms, label
    :rtype: tuple

    """

    calc, label = protocol.get_calculator(
        cpu=cpu, charge=conf.charge, mult=conf.mult, mode="energy"
    )

    atoms = conf.get_ase_atoms(calc)

    if "opt" in protocol.functional.lower():
        with open(f"{label}.{regex_parsing[calc]['ext']}") as f:
            fl = f.readlines()

        geom = get_opt_geometry(fl, protocol.calculator, log)
        set_last_geometry(conf, geom)

    try:
        atoms.get_potential_energy()
    except Exception as e:
        pass

    return atoms, label


def set_last_geometry(conf, geometry):
    """
    Set the last geometry calculated

    :param conf: the conformer instance to optimize
    :type conf: Conformer
    :param geometry: XYZ geometry
    :type geometry: 2D array

    :rtype: None

    """
    conf.last_geometry = geometry[:]
    return None
