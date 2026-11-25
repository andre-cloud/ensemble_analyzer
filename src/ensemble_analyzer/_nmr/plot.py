"""Comparison plots of the NMR spectra of every protocol, one figure per nucleus.

Reads ``NMR_<isotope>_p<N>_comp.xy`` and, when present, the experimental ``nmr_<isotope>_ref.dat``
(two columns: ppm, intensity). Saves ``NMR_<isotope>_comparison.png`` and ``.pickle`` (editable
with ``enan_graph_editor``), with the δ axis decreasing to the right as usual for NMR.
"""

import pickle
import re
from pathlib import Path
from typing import Dict, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

_FILE = re.compile(r"^NMR_(\d+[A-Z][a-z]?)_p(\d+)_comp\.xy$")


def nmr_spectrum_files(idxs: Optional[Sequence[int]] = None, folder: str = ".") -> Dict[str, Dict[int, Path]]:
    out: Dict[str, Dict[int, Path]] = {}
    for f in Path(folder).glob("NMR_*_comp.xy"):
        m = _FILE.match(f.name)
        if not m:
            continue
        nuc, proto = m.group(1), int(m.group(2))
        if idxs is not None and proto not in [int(i) for i in idxs]:
            continue
        out.setdefault(nuc, {})[proto] = f
    return out


def plot_nmr_comparisons(log, idxs: Optional[Sequence[int]] = None, show: bool = False) -> None:
    for nuc, files in sorted(nmr_spectrum_files(idxs).items()):
        fig, ax = plt.subplots(figsize=(10, 5))
        x_min, x_max = np.inf, -np.inf
        for proto, path in sorted(files.items()):
            x, y = np.loadtxt(path, unpack=True)
            ax.plot(x, y / np.max(np.abs(y)), lw=1, alpha=0.8, label=f"Protocol {proto}")
            x_min, x_max = min(x_min, x.min()), max(x_max, x.max())
        exp = Path(f"nmr_{nuc.lower()}_ref.dat")
        if exp.exists():
            xe, ye = np.loadtxt(exp, unpack=True)
            ax.plot(xe, ye / np.max(np.abs(ye)), color="black", lw=1.2, label="Experimental")
        ax.set_xlim(x_max, x_min)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel(rf"$\delta$ ({nuc}) / ppm")
        ax.set_ylabel("Intensity")
        ax.legend(fancybox=True)
        fig.tight_layout()
        fname = f"NMR_{nuc}_comparison"
        with open(f"{fname}.pickle", "wb") as pf:
            pickle.dump(fig, pf)
        fig.savefig(f"{fname}.png", dpi=300)
        if log is not None:
            log.info(f"Saved {fname}.png")
        if show:
            plt.show()
        plt.close(fig)
