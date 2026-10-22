# CLI Tools Reference

In addition to the main `ensemble_analyzer` driver, the package provides a suite of command-line interface (CLI) tools to assist with setup, pre-processing configuration, and post-processing analysis.

## 1. Installation Check (`enan_check`)

A diagnostic tool to verify the installation of python dependencies, the linkage to external QM software (ORCA/Gaussian/NWChem), the optional semi-empirical and ML backends, and the user directories of EnAn.

**Usage:**
```bash
enan_check
```

Checks performed:

- Presence of core Python libraries (`numpy`, `ase`, `scipy`, etc.).
- Availability of the orca executable in the system `PATH`.
- Correct configuration of the `ORCAVERSION` environment variable.
- Availability of Gaussian executables (`g16` or `g09`).
- Availability of the NWChem executable.
- Optional backends: TBLite and the ML potentials (AIMNet, UMA, FairChem, MACE).
- `ENAN_MODELS_DIR`: directory of the ML model weights (default `~/.ensemble_analyzer/models/`).
- `ENAN_TEMPLATES_DIR`: directory of the Protocol Wizard user templates (default `~/.ensemble_analyzer/template/`).

Only missing Python dependencies and ORCA issues make the check fail; everything else is reported as a warning.

---

## 2. Protocol Wizard (`enan_protocol_wizard`)
An interactive Text User Interface (TUI) to generate the `protocol.json` file. It guides the user through the creation of computational steps without needing to manually edit JSON files.

**Usage**
```bash
enan_protocol_wizard
```
Features:

- Three Configuration Levels:
  - Basic: Essential parameters only.
  - Intermediate: Adds constraints and internal coordinate monitoring.
  - Advanced: Full control over energetic thresholds (`thrG`, `thrB`, *etc*.) and pruning settings.
- Fuzzy Search: Quickly find functionals and basis sets from the internal database (for ORCA).

The basic and Intermediate levels allow to less-expert users to create from scratch a protocol file, using the default values of some configuration. 

**Templates**

*New from template* opens a selector with two sections:

- **System Templates**: built-in quick-start protocols (SP, OPT, OPT+FREQ, TS search, TD-DFT, …).
- **User Templates**: your own protocols, read from `~/.ensemble_analyzer/template/` (override with the `ENAN_TEMPLATES_DIR` environment variable).

Each user template is a `.json` file in that folder and is listed under its file name (without extension). The file uses the same format as a `protocol.json` written by the wizard, so the easiest way to create one is to build a protocol, save it, and copy it there:

```bash
mkdir -p ~/.ensemble_analyzer/template
cp protocol.json ~/.ensemble_analyzer/template/my_opt_freq.json
```

A plain JSON list of steps (`[{...}, {...}]`) is also accepted. Files that cannot be read are skipped and shown with a warning in the selector.

---

## 3. Spectra Regrapher (`enan_regraph`)
Regenerates spectral plots (IR, VCD, UV, ECD) without re-running the quantum calculations. This is useful when you want to adjust convolution parameters (FWHM, Shift) or change the plotting units (nm vs eV) after the main run.

Usage:

```bash
enan_regraph [protocol_ids ...] [options]
```

Arguments:
- `protocol_ids`: (Required) List of protocol numbers to regraph (e.g., 0 1 2).

Options:
- `-rb`, `--read-boltz` `<int>`: Use Boltzmann populations from a specific protocol step instead of recalculating them.
- `-no-nm`: Do not save graphs with nanometer (nm) x-axis (saves only default units).
- `-w`, `--weight`: Overlay the weighting function used for comparison on the plot.
- `--disable-color`: Disable colored terminal output.

Example: Recalculate spectra for step 1 using new settings in settings.json:

```bash
enan_regraph 1
````
---

## 4. Average Energy Analyzer (`enan_get_energies`)

Calculates weighted average energies (E, H, G) for the ensemble at a specific temperature. It allows for thermodynamic recalculations (different T or P) and arithmetic operations between protocol steps.

Usage:

```bash
enan_get_energies [options]
```
Core Options:
- `-T`, `--temp` `<float>`: Temperature in Kelvin for recalculation (default: 298.15 K).
- `-o`, `--output` `<file>`: Output log filename.
- `--pressure` `<float>`: Pressure in kPa (default: 101.325).

Thermo Parameters:
- `--cut-off` `<float>`: qRRHO cut-off frequency in cm⁻¹ (default: 100.0).
- `--alpha` `<int>`: qRRHO damping factor (default: 4).
- `--linear`: Treat molecules as linear rigid rotors.

Math Operations:
- `--sub` `<p1> <p2>`: Calculate difference Avg(Protocol 1) - Avg(Protocol 2).
- `--add` `<p1> <p2>`: Calculate sum Avg(Protocol 1) + Avg(Protocol 2).

Example: Compare energies between protocol 2 and 1 at 310K:

```bash
enan_get_energies -T 310 --sub 2 1
```

---

## 5. Graph Editor (`enan_graph_editor`)
A tool to modify the appearance of the generated Matplotlib plots (saved as .pickle files). It supports both an interactive TUI and a batch mode for automated editing.

Usage:

```bash
enan_graph_editor <pickle_file> [options]
```
**Interactive Mode (TUI)**

Simply run without arguments to enter the interactive menu:
```bash
enan_graph_editor IR_comparison.pickle
```
*Allows renaming labels, changing colors, line styles, widths, and transparency.*

**Batch Mode**

Use the --batch flag to apply changes via command line arguments.

Options:
- `--list`, `-l`: List current legend labels and exit.
- `--rename`, `-r` `<old> <new>`: Rename a legend label.
- `--color`, `-c` `<label> <color>` : Change line color.
- `--linestyle`, `-ls` `<label> <style>`: Change line style (e.g., -, --, : ).
- `--format`, `-f`: Output format (png, pdf, svg, pickle).


Example: Rename a protocol and change its color to red, saving as PNG:

```bash
enan_graph_editor IR_comparison.pickle --batch \
            --rename "Protocol 1" "B97-3c" \
            --color "B97-3c" "red" \
            --format png
```