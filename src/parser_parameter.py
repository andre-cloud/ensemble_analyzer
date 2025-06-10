import re
import os
import numpy as np

try:
    from src.regex_parsing import regex_parsing
    from src.rrho import free_gibbs_energy
except ImportError as e:  # pragma: no cover
    print(e)
    from regex_parsing import regex_parsing
    from rrho import free_gibbs_energy

EH_TO_KCAL = 627.5096080305927


def get_param(x, calculator, param):
    """
    Parsing for Rotational Constant. Thought for a map function

    :param x: line of the output
    :type x: str
    :param calculator: name of the calculator used
    :type calculator: str
    :param param: the parameter to be parsed
    :type param: str

    :return: line
    :rtype: str

    """
    if re.search(regex_parsing[calculator][param], x):
        return x


def get_freq(fl: str, calc: str) -> np.ndarray:
    """
    Parsing for frequencies

    :param fl: piece of log file
    :type fl: str
    :param calc: calculator type
    :type calc: str

    :return: frequency [cm-1]
    :rtype: np.array
    """

    fl = (
        "\n".join(
            "".join(fl)
            .split(regex_parsing[calc]["s_freq"])[-1]
            .strip()
            .splitlines()[4:]
        )
        .split(regex_parsing[calc]["e_freq"])[0]
        .strip()
        .splitlines()
    )
    freq = np.array(
        [
            float(i.split()[regex_parsing[calc]["idx_freq"]])
            for i in fl
            if float(i.split()[regex_parsing[calc]["idx_freq"]]) != 0.0
        ]
    )
    return freq


def get_opt_geometry(fl: str, calc: str, log) -> np.ndarray:
    """Fetch the geometry from the calculator output

    :param fl: file read
    :type fl: str
    :param calc: calculator used
    :type calc: str
    :param log: logger instance
    :type log: logger

    :return: the geometry (XYZ) for each atom
    :rtype: np.ndarray
    """

    opt_done = regex_parsing[calc]["opt_done"] in fl
    if not opt_done:
        log.error(f"The optimization did not find a stationary point")

    fl = "".join(fl)
    fl = fl.split(regex_parsing[calc]["opt_done"])[-1]
    fl = (
        fl.split(regex_parsing[calc]["geom_start"])[-1]
        .split(regex_parsing[calc]["break"])[0]
        .strip()
        .splitlines()
    )

    if calc == "orca":
        geom = np.array([i.split()[1:] for i in fl if i], dtype=float)

    return geom


def tranform_float(freq):
    """Transform into a float number a string. Thought for a map function

    :param freq: frequency to be transformed
    :type freq: float

    :return: frequency transformed
    :rtype: str
    """
    return f"{freq:.2f}"


def get_conf_parameters(conf, number: int, p, time, temp: float, log) -> bool:
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

    with open(os.path.join(conf.folder, f"protocol_{number}.out")) as f:
        fl = f.readlines()

    try:
        e = float(
            list(filter(lambda x: get_param(x, p.calculator, "E"), fl))[-1]
            .strip()
            .split()[-1]
        )
    except Exception as e:  # pragma: no cover:
        log.error(e)
        return False

    if 'opt' in [p.functional.lower().split()+p.add_input.lower().split()]:
        conf.last_geometry = get_opt_geometry(fl, p.calculator, log)

    freq = np.array([])
    if p.freq or ("freq" in p.add_input):
        freq = get_freq(fl, p.calculator) * p.freq_fact
        log.info(
            f"{conf.number} has {freq[freq<0].size} imaginary frequency(s): {', '.join(list(map(tranform_float, freq[freq<0])))}"
        )
        if freq.size == 0:
            log.error(("\n".join(fl[-6:])).strip())
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nNo frequency present in the calculation output.\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise IOError("No frequency in the output file")

    try:
        B = np.array(
            list(filter(lambda x: get_param(x, p.calculator, "B"), fl))[-1]
            .strip()
            .split(":")[-1]
            .split(),
            dtype=float,
        )
        b = np.linalg.norm(B)
    except Exception: 
        log.warning('\tB not found')
        b = 1

    try:
        M = np.linalg.norm(
            np.array(
                list(filter(lambda x: get_param(x, p.calculator, "m"), fl))[-1]
                .strip()
                .split(":")[-1]
                .split(),
                dtype=float,
            )
        )
    except Exception:
        log.warning('\tM not found')
        M = 1

    g = ""
    if freq.size > 0:
        g = free_gibbs_energy(
                SCF=e, T=temp, freq=freq, mw=conf.weight_mass, B=B, m=p.mult
            )
    else:
        g_e = conf.energies.get(str(int(number)-1), {}).get("G-E")
        if g_e is not None:
            g = e + g_e

    conf.energies[str(number)] = {
        "E": e * EH_TO_KCAL if e else e,  # Electronic Energy [kcal/mol]
        "G": g * EH_TO_KCAL if g else None,  # Free Gibbs Energy [kcal/mol]
        "B": b if b else 1,  # Rotatory Constant [cm-1]
        "m": M if M else 1,  # dipole momenti [Debye]
        "time": time,  # elapsed time [sec]
        "G-E" : (g-e) if g and e else None,  # G-E [Eh]
    }

    return True


def get_data_for_graph(conformers, protocol, log):
    confs = [i for i in conformers if i.active]
    if protocol.freq or "tddft" in protocol.add_input:
        for i in confs:
            with open(
                os.path.join(os.getcwd(), i.folder, f"protocol_{protocol.number}.out")
            ) as f:
                fl = f.read()

            parse_graph(fl, i, protocol, log)

        return True
    return False


def parse_graph(fl, conf, protocol, log):
    reg = {
        "UV": {
            "start": regex_parsing[protocol.calculator]["s_UV"],
            "end": regex_parsing[protocol.calculator]["break"],
            "x": regex_parsing[protocol.calculator]["idx_en_tddft"],
            "y": regex_parsing[protocol.calculator]["idx_imp_tddft"],
        },
        "ECD": {
            "start": regex_parsing[protocol.calculator]["s_ECD"],
            "end": regex_parsing[protocol.calculator]["break"],
            "x": regex_parsing[protocol.calculator]["idx_en_tddft"],
            "y": regex_parsing[protocol.calculator]["idx_imp_tddft"],
        },
        "IR": {
            "start": regex_parsing[protocol.calculator]["s_IR"],
            "end": regex_parsing[protocol.calculator]["break"],
            "x": regex_parsing[protocol.calculator]["idx_en_ir"],
            "y": regex_parsing[protocol.calculator]["idx_imp_ir"],
        },
        "VCD": {
            "start": regex_parsing[protocol.calculator]["s_VCD"],
            "end": regex_parsing[protocol.calculator]["break"],
            "x": regex_parsing[protocol.calculator]["idx_en_ir"],
            "y": regex_parsing[protocol.calculator]["idx_imp_vcd"],
        },
    }

    conf.energies[str(protocol.number)]["graph"] = {}

    for i in reg.keys():
        if reg[i]["start"] not in fl:
            log.error(
                f"Graph of {conf.number} not found in the output file. Saving empty arrays"
            )
            x = np.array([])
            y = np.array([])
        else:
            tmp = (
                fl.split(reg[i]["start"])[-1]
                .split(reg[i]["end"])[0]
                .strip()
                .splitlines()
            )

            x = np.array([float(j.split()[reg[i]["x"]]) for j in tmp])
            y = np.array([float(j.split()[reg[i]["y"]]) for j in tmp])

            conf.energies[str(protocol.number)]["graph"][i] = {
                "x": x,
                "y": y,
            }


if __name__ == "__main__":  # pragma: no cover:
    # class Conf:
    #     def __init__(self, number, mult, folder):
    #         self.number = number
    #         self.mult = mult
    #         self.folder = folder
    #         self.energies = {}

    # c = Conf('1', 1, 'conf_1')

    # get_conf_parameters(c, 0, 1, 298.15, log)
    # print(c.energies)

    from launch import restart
    from logger import create_log
    import json
    from IOsystem import SerialiseEncoder

    log = create_log("test.out")
    ensemble, protocol, start_from = restart()

    print(protocol)

    for i in protocol:
        print(i)
        if get_data_for_graph(ensemble, i, log):
            print(ensemble[0].energies[str(i.number)]["graph"])

    json.dump(
        {i.number: i.__dict__ for i in ensemble},
        open("checkpoint.json", "w"),
        indent=4,
        cls=SerialiseEncoder,
    )

    # if i.freq or "tddft" in i.add_input:
    #     if i.freq:
    #         Computed(ensemble, False, graph_type='IR', protocol=i.number)
    #         Computed(ensemble, False, graph_type='VCD', protocol=i.number)
    #     if "tddft" in i.add_input:
    #         Computed(ensemble, False, graph_type='UV', protocol=i.number)
    #         Computed(ensemble, False, graph_type='ECD', protocol=i.number)
