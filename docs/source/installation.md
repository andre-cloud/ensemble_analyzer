# Installation and Setup

**Ensemble Analyzer (EnAn)** is a Python-based framework that acts as a high-level driver for quantum mechanical calculations. Therefore, a complete installation involves setting up the Python environment and ensuring the external QM software is correctly linked.

## Prerequisites

Before installing the package, ensure your system meets the following requirements:

* **Operating System**: Linux (recommended), macOS, or Windows.
* **Python**: 3.12, 3.13 or 3.14 (tested). The ML backends need PyTorch wheels for your Python version.
* **QM Engine**: A valid installation of ORCA (recommended), Gaussian and/or NWChem.

---

## 1. Python Environment Setup

We recommend a dedicated environment, created with **Conda** (Anaconda, Miniconda, Miniforge) or with `venv`/`uv`.

```bash
# Create and activate a new environment named 'enan'
conda create -n enan "python>=3.12,<3.15"
conda activate enan
```

### Install the package

```bash
pip install ensemble-analyzer
```

This installs EnAn with its core dependencies (enough to drive ORCA, Gaussian and NWChem).

### Optional backends (extras)

| Extra | Installs | Backend |
|-------|----------|---------|
| `[tblite]` | `tblite` | GFN-xTB (semi-empirical) |
| `[aimnet]` | `aimnet` | AIMNet2 |
| `[uma]` / `[fairchem]` | `fairchem-core` | UMA / FAIRChem models |
| `[mace]` | `mace-torch` | MACE |
| `[skala]` | `skala` | Skala XC functional |
| `[ml]` / `[all]` | all of the above, MACE runtime dependencies included | every ML backend |

```bash
pip install "ensemble-analyzer[ml]"
pip install --no-deps mace-torch
```

**MACE next to UMA/FAIRChem.** `mace-torch` declares `e3nn==0.4.4`, while `fairchem-core` requires `e3nn>=0.5`, so pip refuses to install both in one step. MACE models do run with the newer e3nn: EnAn rebuilds them for the installed version (MACE-OFF23 energies agree within 3·10⁻⁵ eV). Hence `[ml]` installs everything except `mace-torch` itself, which is added with `--no-deps`. With [uv](https://docs.astral.sh/uv/) inside a clone of the repository no extra step is needed (`[tool.uv] override-dependencies`). If you only need MACE, `pip install "ensemble-analyzer[mace]"` is enough.

### Developer install

```bash
git clone https://github.com/andre-cloud/ensemble_analyzer.git
cd ensemble_analyzer
pip install -e ".[dev]"
python -m pytest
```

The package version is derived from the git tags (setuptools-scm): clone the full history, not a shallow copy.

### Conda

All core dependencies are available on conda-forge except `prism_pruner`. A conda-forge package for EnAn is in preparation; until then install the dependencies with conda and EnAn with pip:

```bash
conda create -n enan -c conda-forge "python>=3.12,<3.15" numpy scipy matplotlib-base scikit-learn ase numba tabulate inquirerpy textual sella
conda activate enan
pip install ensemble-analyzer
```

## 2. External QM software

EnAn does not include quantum chemistry codes; it drives them. You must configure at least one calculator. 

### Option A: ORCA (Recommended)
1. **Download & install**: obtain the ORCA binaries directly from the official ORCA forum. 
2. **System Path**: ensure the orca executable is in your system `$PATH` variable. 
```bash
which orca
# Should return a path, e.g., /opt/orca/orca
```
3. **Environment Variable (Crucial)**: EnAn requires the specific ORCA version to handle output parsing correctly. You *must* export the `ORCAVERSION` variable. 
Add the following lines to your shell configuration file (e.g. `~/.bashrc` or `~/.zshrc`): 
```bash
# EnAn specific variable
export ORCAVERSION="6.1.0"
```
Failure to set `ORCAVERSION` will cause the calculator and the parser to crash during the protocol execution.

### Option B: Gaussian
If you have a licensed version of Gaussian installed: 
1. Ensure the standard Gaussian environment variables are loaded (`g09.profile` or `g16.profile`)
2. EnAn will detect the `g16`or `g09`executable automatically from the path. 

## 3. Verification
To verify that EnAn is correctly installed and linked to your QM software (ORCA/Gaussian), simply run the built-in diagnostic tool:
```bash
enan_check
```
This command will check Python dependencies, system paths, and the ORCAVERSION environment variable. The expected output looks as follows: 

```
$ enan_check
Running Ensemble Analyzer Installation Check...


-------------------- 1. Checking Python Dependencies --------------------
[PASS] Library found: numpy
[PASS] Library found: scipy
[PASS] Library found: matplotlib
[PASS] Library found: ase
[PASS] Library found: numba
[PASS] Library found: sklearn
[PASS] Ensemble Analyzer package is installed and importable.

-------------------- 2. Checking ORCA Configuration --------------------
[PASS] ORCA executable found at: /opt/orca/6.1.0/orca
[PASS] ORCAVERSION environment variable is set to: 6.1.0

-------------------- 3. Checking Gaussian Configuration --------------------
[PASS] Gaussian 16 found at: /opt/gaussian/16A03/g16

--------------------------------------------------
SUCCESS: Installation looks correct! You are ready to run EnAn.
```



