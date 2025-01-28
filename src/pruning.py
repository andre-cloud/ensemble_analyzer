from ase.build import minimize_rotation_and_translation
from scipy.constants import R
from tabulate import tabulate
import numpy as np


try:
    from src.ioFile import save_snapshot
    from src.logger import DEBUG, ordinal
except ImportError as e:  # pragma: no cover
    print(e)
    from ioFile import save_snapshot
    from logger import DEBUG, ordinal


def cut_over_thr_max(confs: list, number: str, thrGMAX: float, log) -> None:
    """
    Get conformers over the threshold of the amx energy

    :param confs: whole ensemble list
    :type confs: list
    :param number: number of the protocol to be checked
    :type number: str
    :param thrGMAX: maximum relative energy
    :type thrGMAX: float
    :param log: logger instance
    :type log: logger

    :return: updated ensemble
    :rtype: list
    """

    en = []
    for i in confs:
        e_tmp = -np.inf 
        if i.energies.get(str(number)):
            e_tmp = (
                i.energies[str(number)]["G"]
                if i.energies[str(number)]["G"]
                else i.energies[str(number)]["E"]
            )
        en.append(e_tmp)
    ens = np.array([(i, j) for i, j in zip(confs, en) if i.active])
    ens[:, 1] = ens[:, 1] - min(ens[:, 1])

    log.info(
        f"\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol)"
    )
    for i, en in list(ens):
        if en > thrGMAX:
            i.active = False
            log.info(f"{i.number} - {en:.3f}")
    log.info("\n")


def rmsd(check, ref) -> float:
    r"""
    Compute the Root Mean Squared Deviation (**RMSD**) of two geometries

    .. math::
        RMSD = \sqrt{ \frac {1}{N} ||v_{i} - w_{i}|| }

    :param check: conformer to be compared
    :type check: Conformer
    :param ref: reference
    :type ref: Conformer

    :return: RMSD
    :rtype: float
    """
    ref_pos, check_pos = ref.copy(), check.copy()
    minimize_rotation_and_translation(ref_pos, check_pos)
    return np.sqrt(1 / len(ref.get_positions())) * np.linalg.norm(
        np.array(ref_pos.get_positions()) - np.array(check_pos.get_positions())
    )


def dict_compare(check, conf_ref, deactivate=True) -> dict:  # pragma: no cover
    """
    Create a default dictionary for the comparison

    :rtype: dict
    """
    return {
        "Check": check.number,
        "Ref": conf_ref.number,
        "∆E [kcal/mol]": check.get_energy - conf_ref.get_energy,
        "∆B [e-3 cm-1]": np.abs(check.rotatory - conf_ref.rotatory) * 10**3,
        "∆m [Debye]": np.abs(check.moment - conf_ref.moment),
        "RMSD [Å]": rmsd(check.get_ase_atoms(), conf_ref.get_ase_atoms()),
        "Deactivate": deactivate,
    }


def check(check, conf_ref, protocol, controller: dict) -> bool:
    """
    Control conformer against a reference. Both the following asserts MUST be matched to deactivate the conformer:

    * :math:`B_{check} - B_{reference} < thrB`
    * :math:`G_{check} - G_{reference} < thrG`

    :param check: conformer to be compared
    :type check: Conformer
    :param conf_ref: reference
    :type conf_ref: Conformer
    :param protocol: protocol in order to gain the thresholds
    :type protocol: Protocol
    :param controller: tracker of the deactivated conformers
    :type controller: dict

    :return: true if conf have been deactivated.
    :rtype: bool
    """

    if not conf_ref.active:
        return False

    l = len(controller)

    controller[l] = dict_compare(check, conf_ref, deactivate=False)

    if (
        controller[l]["∆E [kcal/mol]"] < protocol.thrG
        and controller[l]["∆B [e-3 cm-1]"] * 10**-3 < protocol.thrB
    ):
        check.active = False
        check.diactivated_by = conf_ref.number
        controller[l] = dict_compare(check, conf_ref)
        return True

    controller.pop(l)

    return False


def refactor_dict(controller: dict) -> dict:
    """
    Refactor dictionaries in order to print the correct table

    ``{'1': {a: 1, b: 2}, '2' : {a: 3, b: 4}} -> {a : {'1': 1, '2': 3}, b : {'1': 2, '2': 4}}``

    :param controller: dictionary to refactor
    :type controller: dict

    :return: refactored dictionary
    :rtype: dict
    """

    if not controller:
        return {}

    keys = list(controller[0].keys())
    d = {i: [] for i in keys}

    for i in controller:
        for j in d:
            d[j].append(controller[i][j])
    return d


def check_ensemble(confs: list, protocol, log) -> list:
    """
    Check the ensemble:

    1. Over energy threshold
    2. Assert if duplicate conformers with energy and B comparison

    :param confs: whole ensemble list
    :type confs: list
    :param protocol: protocol instance
    :type protocol: Protocol
    :param log: logger instance
    :type log: logger

    :return: ensemble pruned with conformers deactivated
    :rtype: list
    """

    if protocol.graph:
        log.info(
            f"Since graph calculation is detected in this part ({ordinal(int(protocol.number))}), PRUNING NOT EXECUTED"
        )
        return confs

    cut_over_thr_max(confs, protocol.number, protocol.thrGMAX, log)

    if DEBUG:
        save_snapshot(f"after_protocol_{protocol.number}_before_check.xyz", confs, log)

    controller = {}

    for idx, i in enumerate(confs):
        if not i.active:
            continue  # Not check the non active conformers
        for j in range(0, idx):
            if check(i, confs[j], protocol, controller):
                break

    controller = refactor_dict(controller)

    log.info("")
    log.info(tabulate(controller, headers="keys", floatfmt=".3f"))
    log.info("")

    return confs


def calculate_rel_energies(conformers: list, T: float) -> None:
    """
    Relative energy

    :param conformers: whole ensemble list
    :type conformers: list
    :param T: temperature [K]
    :type T: float

    :return: None
    """

    c = [i for i in conformers if i.active]
    ens = np.array([i.get_energy for i in conformers if i.active])
    ens -= min(ens)
    bolz = np.exp((-ens * 4186) / (R * T))
    pop = (bolz / np.sum(bolz)) * 100
    for idx, i in enumerate(list(ens)):
        c[idx]._last_energy["Erel"] = i
        c[idx]._last_energy["Pop"] = pop[idx]

    return None
