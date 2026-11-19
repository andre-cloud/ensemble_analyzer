"""Tasks (FAIRChem ``task_name``) and heads (MACE ``head``) supported by known
ML-potential checkpoints, as ``{engine: {task: [model files]}}``.

Sources:
  * Task/head lists read from the checkpoints themselves on 2026-11-18
    (fairchem ``load_predict_unit(...).dataset_to_tasks``,
    mace ``MACECalculator(...).available_heads``).
  * FAIRChem/UMA: https://fair-chem.github.io/core/uma.html (``charge`` and
    ``spin`` are used by ``omol`` only) and https://pypi.org/project/fairchem-core/
    (``oc22``/``oc25`` exist in the 1p2 models only).
  * MACE: an unknown ``head`` silently falls back to the last head
    (https://github.com/ACEsuit/mace/blob/v0.3.17/mace/calculators/mace.py, lines 339-343),
    hence the explicit check here.

Single-head MACE models expose a head called ``default`` or ``Default``
(depending on the release); heads are matched case-insensitively.
Checkpoints not listed here are accepted with a "no sanity check" warning.
"""

import warnings
from pathlib import Path

_UMA_ALL = ["uma-s-1.pt", "uma-s-1p1.pt", "uma-s-1p2.pt", "uma-m-1p1.pt"]

_EGRET = ["EGRET_1.model", "EGRET_1E.model", "EGRET_1M.model", "EGRET_1S.model", "EGRET_1T.model"]
_MACE_MH = ["mace-mh-0.model", "mace-mh-1.model"]

MODEL_TASKS: dict[str, dict[str, list[str]]] = {
    "fairchem": {
        "omol": _UMA_ALL + [
            "uma-s-1p2p1-ft-d4rev.pt",
            "esen_sm_conserving_all.pt", "esen_sm_direct_all.pt", "esen_md_direct_all.pt",
            "AllScAIP-OMol102M-md-cons.pt", "AllScAIP-OMol102M-md-d.pt",
        ],
        "omat": _UMA_ALL + ["uma-s-1p2p1-ft-r2scan.pt"],
        "omc": list(_UMA_ALL),
        "odac": list(_UMA_ALL),
        "oc20": list(_UMA_ALL),
        "oc22": ["uma-s-1p2.pt"],
        "oc25": ["uma-s-1p2.pt"],
    },
    "mace": {
        "omol": ["MACE-omol-0-extra-large-1024.model"] + _MACE_MH,
        "omat_pbe": ["mace-omat-0-small.model"] + _MACE_MH,
        "matpes_r2scan": list(_MACE_MH),
        "mp_pbe_refit_add": list(_MACE_MH),
        "spice_wB97M": list(_MACE_MH),
        "oc20_usemppbe": list(_MACE_MH),
        "rgd1_b3lyp": ["mace-mh-0.model"],
        "default": _EGRET + [
            "mace-large-density-agnesi-stress.model",
            "mace-medium-density-agnesi-stress.model",
            "mace-small-density-agnesi-stress.model",
            "mace-mp-0b3-medium.model", "mace-mpa-0-medium.model", "mace-omat-0-medium.model",
            "MACE-matpes-pbe-omat-ft.model", "MACE-matpes-r2scan-omat-ft.model",
        ],
        "Default": [
            "MACE-OFF23_small.model", "MACE-OFF23_medium.model", "MACE-OFF23_large.model",
            "MACE-POLAR-1-S.model", "MACE-POLAR-1-M.model", "MACE-POLAR-1-L.model",
            "MACE_MPtrj_2022.9.model", "mace_agnesi_medium.model", "mace_agnesi_small.model",
            "2023-12-03-mace-128-L1_epoch-199.model", "2023-12-10-mace-128-L0_energy_epoch-249.model",
        ],
    },
}

# "uma" is the same fairchem-core backend
MODEL_TASKS["uma"] = MODEL_TASKS["fairchem"]


def tasks_for_model(engine: str, model: str) -> list[str] | None:
    """Tasks/heads of a known checkpoint, or None if it is not listed.

    Matches the file name; a prefixed copy (e.g. the mace download cache name
    ``mace-polar-polar-1-m-MACE-POLAR-1-M.model``) matches by suffix.
    """
    name = Path(str(model)).name
    found = [
        task for task, models in MODEL_TASKS.get(engine, {}).items()
        if any(name == m or name.endswith("-" + m) for m in models)
    ]
    return found or None


def check_model_task(engine: str, model: str, task: str | None) -> None:
    """Raise if ``task`` is not supported by a known checkpoint; warn if unknown.

    ``task=None`` means "model default" and is left to the backend, which
    already errors when a multi-head model has no default head.
    """
    tasks = tasks_for_model(engine, model)
    if tasks is None:
        warnings.warn(
            f"{engine} model '{Path(str(model)).name}' is not in the EnAn model table: "
            f"task={task!r} accepted with no sanity check.",
            stacklevel=3,
        )
        return
    if task is None:
        return
    if engine == "mace":
        ok = task.lower() in {t.lower() for t in tasks}
    else:
        ok = task in tasks
    if not ok:
        raise ValueError(
            f"{engine} model '{Path(str(model)).name}' does not support task {task!r}. "
            f"Supported: {', '.join(sorted(tasks))}."
        )
