from tabulate import tabulate

import numpy as np
from src.constants import *

try:
    from src.ioFile import save_snapshot
    from src.logger import DEBUG, ordinal
    from src.conformer import Conformer
except ImportError as e:  # pragma: no cover
    print(e)
    from ioFile import save_snapshot
    from logger import DEBUG, ordinal
    from conformer import Conformer


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

    ens = np.array([(j, (j.energies[str(number)]["G"]
                if not np.isnan(j.energies[str(number)]["G"])
                else j.energies[str(number)]["E"])) for j in confs if i.active])
    
    ens[:, 1] = ens[:, 1] - min(ens[:, 1])

    log.info(
        f"\nGetting number of conformers lying out of the energy windows (over {thrGMAX} kcal/mol)"
    )
    for i, en in list(ens):
        if en > thrGMAX:
            i.active = False
            log.info(f"{i.number} - {en:.3f}")
    log.info("\n")


def rmsd(check, ref, include_H=False) -> float:
    r"""
    Compute the Root Mean Squared Deviation (**RMSD**) of two geometries based on the distance matrix's eigenvalues.

    .. math::
        RMSD = \sqrt{ \frac {1}{N} \left(\lambda_{i} - lambda_{i}\right)^2 }

    :param check: conformer to be compared
    :type check: Conformer
    :param ref: reference
    :type ref: Conformer

    :return: RMSD
    :rtype: float
    """
    ref_pos, check_pos = ref.distance_matrix(include_H), check.distance_matrix(
        include_H
    )
    eva_ref, _ = np.linalg.eig(ref_pos)
    eva_check, _ = np.linalg.eig(check_pos)
    return np.sqrt(1 / len(eva_ref) * np.sum((eva_ref - eva_check) ** 2))


def dict_compare(
    check, conf_ref, include_H, deactivate=True, RMSD=None
) -> dict:  # pragma: no cover
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
        "λi RMSD": RMSD if RMSD else rmsd(check, conf_ref, include_H),
        "Deactivate": deactivate,
    }


def check_dE_dB(check, conf_ref, protocol, controller: dict, include_H: bool) -> bool:
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
    if check.rotatory == 1:
        return False

    l = len(controller)

    controller[l] = dict_compare(check, conf_ref, deactivate=False, include_H=include_H)

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


def check_enantiomer(check, ref, protocol, controller: dict) -> bool:
    """Check if two conformers are enatiomeric conformations

    Args:
        :param check: conformer to be compared
        :type check: Conformer
        :param ref: reference
        :type ref: Conformer
        :param protocol: protocol in order to gain the thresholds
        :type protocol: Protocol
        :param controller: tracker of the deactivated conformers
        :type controller: dict

    Returns:
        :return: Returns True if the two conformers are enantiomers
        :type: bool
    """

    l = len(controller)

    RMSD = rmsd(check, ref, include_H=False)
    controller[l] = dict_compare(
        check, ref, deactivate=False, include_H=False, RMSD=RMSD
    )

    if RMSD < protocol.thrRMSD_enantio:
        check.active = False
        check.diactivated_by = ref.number
        controller[l] = dict_compare(check, ref)
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


def check_ensemble(
    confs: list, protocol, log, include_H
) -> list:  # exclude_enantiomers
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
            f"Graph calculation is detected in this part ({ordinal(int(protocol.number))}), PRUNING NOT EXECUTED"
        )
        return confs

    if protocol.no_prune:
        log.info(
            f"No prune setting is detected in this part ({ordinal(int(protocol.number))}), PRUNING NOT EXECUTED"
        )
        return confs

    # Energy Thresholds
    cut_over_thr_max(confs, protocol.number, protocol.thrGMAX, log)

    if DEBUG:
        save_snapshot(f"after_protocol_{protocol.number}_before_check.xyz", confs, log)

    controller = {}

    # OneToOne CHECK
    for idx, check in enumerate(confs):
        if not check.active:
            continue  # Not check the non active conformers

        mirror = check.last_geometry.copy()
        mirror[..., 2] *= -1
        check_mirror = Conformer.load_raw(check.__dict__)
        check_mirror.last_geometry = mirror

        for j in range(0, idx):

            # if exclude_enantiomers:
            #     if check_enantiomer(check=check_mirror, conf_ref=confs[j], controller=controller):
            #         break

            if check_dE_dB(
                check=check,
                conf_ref=confs[j],
                protocol=protocol,
                controller=controller,
                include_H=include_H,
            ):
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
