# Protocol Configuration Guide

The protocol file is the core configuration for **Ensemble Analyzer**. It is structured as a JSON dictionary where each key (e.g., `"0"`, `"1"`) represents a sequential computational step.

## Core Computational Keywords

These parameters define the level of theory and the type of calculation to be performed by the QM engine.

* **`functional`** (str, *Required*)
    The DFT functional or semi-empirical method to use (e.g., `"B97-3c"`, `"wB97X-D4"`, `"xtb"`).

* **`basis`** (str)
    The basis set definition. If using composite methods (like `r2SCAN-3c`), this is automatically handled or can be omitted.

* **`opt`** (bool)
    If `true`, performs a geometry optimization for the current step.

* **`freq`** (bool)
    If `true`, performs a frequency calculation. This enables vibrational analysis and qRRHO thermochemical corrections.

* **`charge`** (int)
    The total charge of the system (default: `0`).

* **`mult`** (int)
    The spin multiplicity of the system (default: `1`).

* **`solvent`** (dict)
    Configuration for implicit solvation models.
    * `solvent` (str): Name of the solvent (e.g., `"water"`, `"chcl3"`).
    * `smd` (bool): If `true`, uses the **SMD** solvation model; otherwise, uses **CPCM**.

## Advanced Calculation Control

* **`add_input`** (str)
    Additional keywords or blocks passed directly to the external QM engine (ORCA/Gaussian) input file. **No sanity check performed**. 

* **`read_orbitals`** (int)
    Specifies the index of a previous step to read orbitals/guess from (e.g., `"0"`). Useful for SCF convergence in difficult cases.

* **`skip_opt_fail`** (bool)
    If `true`, conformers that fail to converge during optimization are automatically deactivated instead of crashing the workflow.

* **`monitor_internals`** (list)
    A list of atom indices to track specific internal coordinates in the log output.
    * *Example:* `[[0, 1], [2, 3, 4]]` monitors a bond length and an angle.

* **`block_on_retention_rate`** (bool)
    If `true`, the program will halt execution if the number of surviving conformers drops below a safety threshold (default 20%), preventing total loss of the ensemble.

## Transition State Analysis

* **`ts`** (bool)
    If `true`, enables transition state optimization. QM calculators add TS-specific keywords (`OptTS` for ORCA, `opt=(ts,calcfc,noeigentest)` for Gaussian); ML calculators use the Sella optimizer instead of LBFGS. Implies `"opt": true`.

* **`loc_freq`** (list of list of int)
    List of atom index groups that the imaginary frequency should localize on. When set, gates the B.1–B.6 post-processing logic: validates whether each significant negative mode falls on these groups (auto-named `Frag1`, `Frag2`, …), displaces and re-optimizes when a spurious mode is detected. Without `loc_freq`, even with `"ts": true`, the TS-specific post-processing is skipped.
    * *Example:* `[[0, 1, 2, 3], [4, 5, 6]]` defines two fragments.

* **`min_overlap`** (float)
    Minimum percentage of squared atomic displacement that must localize on the `loc_freq` fragments for a negative frequency to be considered "on target". Default: `50.0`.

* **`auto_displace`** (bool)
    If `true`, automatically displaces the geometry along the dominant imaginary mode and re-optimizes. Applies to both TS and non-TS optimizations with imaginary frequencies.

* **`displace_scale`** (float)
    Displacement scale factor in Ångström. Default: `0.3`.

* **`neg_freq_threshold`** (float)
    Imaginary frequencies with |ν| ≤ this value (in cm⁻¹) are classified as noise and ignored. Default: `20.0`.

## Refinement & Pruning Settings

* **`cluster`** (int | bool)
    Controls the unsupervised clustering of conformers.
    * If an **integer > 1**: Performs K-Means clustering to reduce the ensemble to that exact number of structures.
    * If **`true`**: Performs clustering with automatic detection of the optimal number of clusters ($k$).

* **`no_prune`** (bool)
    If `true`, completely disables energy and geometric pruning for this specific step.

* **`thrG`** / **`thrB`** (float)
    Overrides the default thresholds for identifying duplicates:
    * `thrG`: Maximum energy difference ($\Delta E$) [kcal/mol].
    * `thrB`: Maximum difference in Rotational Constants ($\Delta B$) [cm⁻¹].

* **`thrGMAX`** (float)
    Overrides the maximum energy window cut-off. Conformers with $\Delta E > \text{thrGMAX}$ (relative to the global minimum) are discarded.