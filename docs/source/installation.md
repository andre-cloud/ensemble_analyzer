# Installation and Setup

**Ensemble Analyzer (EnAn)** is a Python-based framework that acts as a high-level driver for quantum mechanical calculations. Therefore, a complete installation involves setting up the Python environment and ensuring the external QM software is correctly linked.

## Prerequisites

Before installing the package, ensure your system meets the following requirements:

* **Operating System**: Linux (recommended), macOS, or Windows.
* **Python**: Version 3.11 or higher.
* **QM Engine**: A valid installation of ORCA (recommended) and/or Gaussian.

---

## 1. Python Environment Setup

We strongly recommend using **Conda** (via Anaconda or Miniconda) to manage dependencies and avoid conflicts with system libraries.

### Create a Virtual Environment

```bash
# Create a new environment named 'enan'
conda create -n enan python>=3.11

# Activate the environment
conda activate enan
```
### Install the package
You can install EnAn directly using `pip`. This automatically resolve core dependencies.
```bash
pip install ensemble-analyzer
```

For developers who want to modify/add some feature
```bash
git clone
cd ensemble_analyzer
pip install -e .
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



