# Ensemble Analyzer
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Conformer Ensemble Pruning Software


![logo](logo.png)

**EnAn** (Ensemble Analysis) is a Python framework for automated conformational analysis and ensemble processing in computational chemistry workflows.

---

## 🎯 Key Features

### Core Capabilities
- ⚡ **Multi-Protocol Workflows**: Sequential optimization/frequency calculations with automatic pruning
- 🔬 **Quantum Chemistry Integration**: Support for ORCA and Gaussian
- 📊 **Advanced Clustering**: PCA-based conformer clustering with multiple feature extraction methods
- 🎨 **Spectral Analysis**: Generate weighted IR, VCD, UV-vis, and ECD spectra
- 🔄 **Checkpoint System**: Automatic restart capability with atomic file operations

### Analysis Tools
- **Thermochemistry**: Grimme's qRRHO implementation for every integration
- **Clustering Method**: Distance matrix eigenvalues (rotation/translation invariant).
- **Pruning**: Intelligent energy-based filtering with Boltzmann weighting


---

## 📦 Installation

### Requirements
- **Python** ≥ 3.9
- **Core**: NumPy, SciPy, Matplotlib, ASE
- **Clustering**: scikit-learn
- **Acceleration**: Numba

```bash
pip install ensemble-analyzer
```

- Install [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) from the ORCA Forum
- *Optional*: Install Gaussian, if licensed

- Export ORCA verion
```bash
export ORCAVERSION="x.y.z"
```

## 🚀 Quick Start

### 1. Basic Usage
```bash
ensemble_analyzer --ensemble conformers.xyz --protocol protocol.json --output calculation.out --cpu 8 --temperature 298.15
```

### 2. Restart from Checkpoint
```bash
# Automatically resumes from last completed protocol
ensemble_analyzer --restart
```

### 3. Define your protocol file
Create `protocol.json`
```json
{
    "0": {"funcional": "r2SCAN-3c", "opt": true, "freq": true,"cluster": 5, "comment": "Initial Optimization cluster into 5 families"},
    "1": {"funcional": "wB97X-D4rev", "basis": "def2-QZVPPD", "comment": "Single Point energy evaluation"}
}

```
---

### Protocol Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| **Calculation Settings** ||||
| `functional` | str | DFT functional or method | `"B3LYP"`, `"xtb"`, `"HF-3c"` |
| `basis` | str | Basis set (auto for composite methods) | `"def2-SVP"`, `"def2-TZVP"` |
| `calculator` | str | QM program | `"orca"` (default), `"gaussian"` |
| `opt` | bool | Optimize geometry | `true`, `false` |
| `freq` | bool | Calculate frequencies | `true`, `false` |
| `mult` | int | Spin multiplicity | `1` (singlet), `2` (doublet) |
| `charge` | int | Molecular charge | |
| **Pruning Thresholds** ||||
| `thrG` | float | Energy similarity threshold [kcal/mol] | `3.0`, `5.0` |
| `thrB` | float | Rotatory constant threshold [cm⁻¹] | `30.0`, `50.0` |
| `thrGMAX` | float | Energy window cutoff [kcal/mol] | `10.0` |
| `cluster` | bool/int | Enable clustering | `true` (auto), `5` (fixed) |
| `no_prune` | bool | Disable pruning | `false` (default) |
| **Advanced** ||||
| `solvent` | dict | Implicit solvation | `{"solvent": "water", "model": "SMD"}` |
| `constrains` | list | Geometry constraints (only on cartesians)| `[1,2]` (fix cartesians) |
| `monitor_internals` | list | Track bond/angle/dihedral | `[[0,1], [0,1,2]]` |
| `skip_opt_fail` | bool | Skip failed optimizations | `false` (default) |

---
## 🤝 Contributing

### Development Workflow
```bash
# Fork and clone
git clone https://github.com/your-username/enan.git
cd enan

# Create feature branch
git checkout -b feature/awesome-feature

# Install dev dependencies
pip install -e .

# Commit and push
git commit -m "Add awesome feature"
git push origin feature/awesome-feature
```
---

## 📄 License

This software is licensed under the MIT-3.0 License. See the LICENSE file for details.

## 📞 Support and contact

For any questions or support, please contact [by email](mailto:andrea.pellegrini15@unibo.it,paolo.righi@unibo.it).
