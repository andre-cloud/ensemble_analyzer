# CLI Tools Reference

The main driver `ensemble_analyzer` accepts various command-line arguments to configure the runtime environment, physical constants, and post-processing options.

## Input & Flow Control

* **`-e`**, **`--ensemble`** (path)
    The input XYZ file containing the initial conformational ensemble.

* **`-p`**, **`--protocol`** (path)
    Path to the `protocol.json` file. If omitted, the internal default protocol is used.

* **`--restart`**
    Resumes the calculation from the last completed protocol step using the generated `checkpoint.json`.

* **`-o`**, **`--output`** (str)
    Name of the main log file (default: `output.out`).

### Physics & Thermodynamics

* **`-T`**, **`--temperature`** (float)
    Temperature in Kelvin used for Boltzmann weighting and thermochemical analysis (default: `298.15` K).

* **`-P`**, **`--pressure`** (float)
    Pressure in kPa (default: `101.325`).

* **`--linear`**
    Flag to indicate that the molecules should be treated as linear rigid rotors (affects entropy calculations).

* **`-cut-off`**, **`--cut-off`** (float)
    The frequency cut-off value ($\omega_0$) for the Grimme qRRHO approximation (default: `100.0` cm⁻¹).

* **`-no-H`**, **`--exclude-H`**
    If set, excludes Hydrogen atoms from RMSD and PCA distance matrix calculations.


* **`-cpu`** (int)
    Number of CPU cores to allocate for each QM calculation step (default: `1`).

* **`-calc`**, **`--calculator`** (str)
    Selects the external QM engine driver. Options: `"orca"` (default), `"gaussian"`.

## Spectral Analysis (Graphing)

* **`--fwhm-vibro`** / **`--fwhm-electro`**
    Define the Full Width at Half Maximum (FWHM) for convolution.
    * Can be a single float (fixed value).
    * Can be a range `[min, max]` for auto-optimization against a reference spectrum.

* **`--shift-vibro`** / **`--shift-electro`**
    Define the spectral shift parameter.
    * **Vibronic (IR/VCD)**: Multiplicative scaling factor.
    * **Electronic (UV/ECD)**: Additive shift [eV].

* **`--invert`**
    Inverts the y-axis for chiral spectra (VCD/ECD).

* **`--definition`** (int)
    Resolution of the spectral grid, defined as $10^N$ points (default: `4`, i.e., 10,000 points).

## Utilities

* **`-h`**, **`--help`**
    Shows the main help message.

* **`-h-p`**, **`--help-protocol`**
    Prints a JSON template for the protocol file.

* **`-h-t`**, **`--help-threshold`**
    Prints a JSON template for the threshold file.

* **`--disable-color`**
    Disables ANSI colored output in the terminal logs.