from __future__ import annotations

from typing import Any

DEFAULTS: dict[str, Any] = {
    "calculator": "orca",
    "mult": 1,
    "charge": 0,
    "freq_fact": 1.0,
    "maxstep": 0.2,
    "no_prune": False,
    "cluster": False,
    "opt": False,
    "freq": False,
    "ts": False,
    "tda": False,
    "skip_opt_fail": False,
    "block_on_retention_rate": False,
    "auto_displace": False,
    "min_localization": 40.0,
    "displace_scale": 0.3,
    "neg_freq_threshold": 20.0,
    "fmax": 0.01,
    "maxiter": 100000000,
    "nroots": None,
    "thrG": None,
    "thrB": None,
    "thrGMAX": None,
}

ML_CALCULATORS = {"tblite", "aimnet", "mace", "uma", "fairchem", "skala"}

CALCULATOR_CHOICES: list[tuple[str, str]] = [
    ("ORCA", "orca"),
    ("Gaussian", "gaussian"),
    ("NWChem", "nwchem"),
    ("tblite (GFN-xTB)", "tblite"),
    ("ANI/aimnet", "aimnet"),
    ("MACE-MP", "mace"),
    ("UMA (FairChem)", "uma"),
    ("FairChem (OC-20)", "fairchem"),
    ("Skala (Microsoft)", "skala"),
]
