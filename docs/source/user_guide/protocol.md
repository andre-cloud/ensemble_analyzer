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

* **`calculator`** (str)
    Selects the external QM engine driver or ML potential. Options: `"orca"` (default), `"gaussian"`, `"nwchem"`, `"tblite"`, `"aimnet"`, `"uma"`, `"fairchem"`, `"mace"`.

* **`task`** (str)
    Selects the level of theory inside a multi-task ML model. Leave it out to use the model default. Its meaning depends on the `calculator`:
    * `"uma"`, `"fairchem"`: the FAIRChem `task_name` (e.g. `"omol"`, `"omat"`, `"oc20"`, `"odac"`, `"omc"`). Default: `"omol"`, the only task meant for isolated molecules (it uses `charge` and `mult`). The other tasks are trained on periodic systems with different DFT references, so their energies are not comparable with `"omol"`. Single-task checkpoints must be paired with their own task (e.g. `uma-s-1p2p1-ft-r2scan.pt` only supports `"omat"`); an unsupported task raises an error.
    * `"mace"`: the `head` of a multi-head MACE model. Default: the model's default head. Single-head models (e.g. MACE-MP, MACE-OFF) need no `task`.
    * `"aimnet"`: not supported. AIMNet2 models have a single level of theory, so setting `task` raises an error when the protocol is loaded; pick the model via `functional` instead.
    * QM engines (`"orca"`, `"gaussian"`, `"nwchem"`, ...): ignored.

    ```json
    "1": {"calculator": "uma", "functional": "uma-s-1p2p1.pt", "task": "omol", "opt": true}
    ```

* **`add_input`** (str)
    Additional keywords or blocks passed directly to the external QM engine (ORCA/Gaussian) input file. **No sanity check performed**. 

* **`read_orbitals`** (int)
    Specifies the index of a previous step to read orbitals/guess from (e.g., `"0"`). Useful for SCF convergence in difficult cases.

* **`read_population`** (str)
    Specifies the index of a previous step to read Boltzmann populations from.

* **`skip_opt_fail`** (bool)
    If `true`, conformers that fail to converge during optimization are automatically deactivated instead of crashing the workflow.

* **`monitor_internals`** (list)
    A list of atom indices to track specific internal coordinates in the log output.
    * *Example:* `[[0, 1], [2, 3, 4]]` monitors a bond length and an angle.

* **`constraints`** (list[list[int]])
    A list of atom indices groups to constrain during geometry optimization.
    * *Example:* `[[1, 2], [3, 4, 5]]` restrains bond between 1 and 2, and angle between 3, 4, 5.

* **`validators`** (list[list])
    A list of custom regex patterns to validate the output file. Each entry should be `[pattern, expected_value, threshold]`.

* **`comment`** (str)
    A custom comment or description for this protocol step.

* **`block_on_retention_rate`** (bool)
    If `true`, the program will halt execution if the number of surviving conformers drops below a safety threshold (default 20%), preventing total loss of the ensemble.

* **`freq_fact`** (float)
    Scale factor applied to calculated frequencies.

* **`fmax`** (float)
    Convergence threshold for ML optimizers (BFGS/Sella) [eV/Å] (default: `0.01`).

* **`maxstep`** (float)
    Maximum step size for ML optimizers [Å] (default: `0.2`).

* **`maxiter`** (int)
    Maximum number of optimization iterations.

## TD-DFT Settings (UV-Vis & ECD Spectra)

To compute excited states and generate electronic spectra (UV-Vis and ECD), you must define the `nroots` keyword in the protocol step. 

> [!WARNING]
> **Limitation:** TD-DFT calculations are currently **only supported by QM engines** (ORCA, Gaussian, NWChem). They are **NOT supported** when using ML Potentials (MLIPs like AIMNet, MACE) or semi-empirical methods like `tblite`.

* **`nroots`** (int)
    Number of excited states (roots) to calculate. Adding this keyword automatically triggers the TD-DFT module.
* **`tda`** (bool)
    Toggle the Tamm-Dancoff approximation. By default, ORCA uses `true` (TDA on). Set to `false` for full TD-DFT.

### Example: TD-DFT Protocol Step
This example shows a typical workflow: a ground-state optimization and frequency calculation, followed by a single-point TD-DFT calculation (computing 30 excited states).

```json
{
    "0": {
        "calculator": "orca",
        "functional": "wB97X-D4rev",
        "basis": "def2-SVP",
        "opt": true,
        "freq": true
    },
    "1": {
        "calculator": "orca",
        "functional": "wB97X-D4rev",
        "basis": "def2-TZVPPD",
        "nroots": 30,
        "tda": false,
        "comment": "TD-DFT single point for UV-Vis/ECD generation"
    }
}
```

## Transition State Analysis

* **`ts`** (bool)
    If `true`, enables transition state optimization. QM calculators add TS-specific keywords (`OptTS` for ORCA, `opt=(ts,calcfc,noeigentest)` for Gaussian); ML calculators use the Sella optimizer instead of LBFGS. Implies `"opt": true`.

* **`loc_freq`** (list of list of int)
    List of atom index groups that the imaginary frequency should localize on. When set, gates the B.1–B.6 post-processing logic: validates whether each significant negative mode falls on these groups (auto-named `Frag1`, `Frag2`, …), displaces and re-optimizes when a spurious mode is detected. Without `loc_freq`, even with `"ts": true`, the TS-specific post-processing is skipped.
    * *Example:* `[[0, 1, 2, 3], [4, 5, 6]]` defines two fragments.

* **`min_localization`** (float)
    Minimum percentage of squared atomic displacement that must localize on the `loc_freq` fragments for a negative frequency to be considered "on target". Default: `50.0`.

* **`auto_displace`** (bool)
    If `true`, automatically displaces the geometry along the dominant imaginary mode and re-optimizes. Applies to both TS and non-TS optimizations with imaginary frequencies.

* **`displace_scale`** (float)
    Displacement scale factor in Ångström. Default: `0.3`.

* **`neg_freq_threshold`** (float)
    Imaginary frequencies with |ν| ≤ this value (in cm⁻¹) are classified as noise and ignored. Default: `20.0`.

### Example: Transition State Optimization
This step optimizes a transition state using Sella (if using an ML calculator) or standard OptTS (ORCA), and automatically displaces/re-optimizes the geometry if the transition state frequency doesn't involve the specified breaking/forming bonds between atoms 1-4 and 5-10.

```json
{
    "0": {
        "calculator": "orca",
        "functional": "wB97X-D4rev",
        "basis": "def2-TZVP",
        "ts": true,
        "loc_freq": [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9]],
        "auto_displace": true,
        "min_localization": 50.0
    }
}
```

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

* **`prune_method`** (str, default `"rotconst"`)
    Selects the duplicate-removal backend (the `thrGMAX` energy window is always applied first):
    * `"rotconst"`: energy + rotational-constant norm (`thrG`, `thrB`). Fast, no extra dependency.
    * `"prism"`: [PRISM Pruner](https://github.com/ntampellini/prism_pruner) — heavy-atom Kabsch RMSD + maximum-deviation filter (optionally preceded by a moment-of-inertia filter, see `prism_moi`).
    * `"prism_rotcorr"`: as `"prism"`, plus an RMSD corrected for degenerate rotamers (e.g. methyl / *t*-Bu rotations). Slowest.

    > **Note on `rotconst` for large molecules.** `thrB` is an absolute threshold, while rotational constants shrink with molecular size. In a test on raw CREST ensembles with the default `opt` thresholds (`thrG` 0.5 kcal/mol, `thrB` 5·10⁻⁵ cm⁻¹), `rotconst` removed 541/673 conformers of a 71-atom system and 651/675 of a 220-atom system; most removed conformers differed from their reference by > 0.5 Å heavy-atom RMSD, also after rotamer correction. For large or flexible systems prefer `prism` or tighten `thrB`.

    With PRISM, only pairs with $|\Delta E| <$ `thrG` are compared; `thrB` is not used. The lowest-energy member of each duplicate pair is always retained. `prism_pruner` is installed with EnAn as a core dependency.

* **`prism_moi`** (bool, default `false`)
    Run PRISM's moment-of-inertia stage before the RMSD stage. It is fast and reflection-invariant (mirror-image conformers are merged, as with `rotconst`), but with the PRISM default `prism_moi_dev = 0.01` it merged clearly distinct conformers of large molecules in our tests (71- and 220-atom ensembles), whereas the RMSD stage alone produced no such merges. If enabled, a tighter `prism_moi_dev` (e.g. `0.001`) is recommended.

* **`prism_moi_dev`** (float, default `0.01`)
    Relative deviation threshold on all three principal moments (0.01 = 1 %).

* **`prism_max_rmsd`** / **`prism_max_dev`** (float, defaults `0.25` / `2 × prism_max_rmsd`) [Å]
    Two conformers are duplicates if heavy-atom RMSD < `prism_max_rmsd` **and** the maximum atomic deviation < `prism_max_dev`. Hydrogens are included only if the run uses `include_H`.

* **`prism_timeout`** (int, default `60`) [s]
    Timeout of each PRISM stage; on timeout the stage returns what it has pruned so far.

```json
{
    "0": {"calculator": "tblite", "functional": "GFN2-xTB", "opt": true,
          "prune_method": "prism", "comment": "xTB pre-optimisation, PRISM de-duplication"},
    "1": {"functional": "r2SCAN-3c", "opt": true, "freq": true,
          "prune_method": "prism_rotcorr", "prism_max_rmsd": 0.2,
          "prism_moi": true, "prism_moi_dev": 0.001}
}
```