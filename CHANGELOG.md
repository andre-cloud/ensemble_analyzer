# Changelog

All notable changes to this project will be documented in this file.


## [Unreleased]

### Added
- **Transition State (TS) Optimizer & Mode Validation**: Automated fragment-based localization of imaginary frequencies to validate transition states. Includes automatic displacement along spurious imaginary modes and subsequent re-optimization.
- **TD-DFT Support**: Added full support for time-dependent density functional theory (TD-DFT) calculations, including Tamm-Dancoff Approximation (TDA) control (`tda`) and configuration of the number of excited roots (`nroots`).
- **Electronic Spectral Generation**: Implemented Gaussian line-shape convolution for UV-Vis and Electronic Circular Dichroism (ECD) spectra generation.
- CLI Graph Editor (`enan_graph_editor`) and Spectra Regrapher (`enan_regraph`) tools.
- **ML Potentials Support**: Integration with AIMNet, UMA, FairChem, and MACE via the ASE calculator wrapper.
- **`task` protocol keyword**: selects the level of theory of multi-task ML models — FAIRChem `task_name` for `uma`/`fairchem` (default `omol`), `head` for multi-head `mace` models (default: model head). Setting it with `aimnet` raises an error at protocol load. Also exposed in the Protocol Wizard.
- Automatic Boltzmann weighting and conformational pruning strategies based on Euclidean Distance Matrix eigenvalues.
- **User templates in the Protocol Wizard**: protocols saved as `.json` in `~/.ensemble_analyzer/template/` (or `ENAN_TEMPLATES_DIR`) appear under *User Templates* next to the built-in *System Templates*. The folder layout is free (sub-folders and symlinked group folders are shown as sub-sections), and the new *Save as template* button stores the current protocol there.

- **PRISM Pruner backend** (`prism_pruner` is a core dependency): new protocol keyword `prune_method` (`rotconst` default | `prism` | `prism_rotcorr`) with `prism_max_rmsd`, `prism_max_dev`, `prism_moi`, `prism_moi_dev`, `prism_timeout`. Exposed in the Protocol Wizard and checked by `enan_check`. EnAn pre-sorts by energy (works around a re-sorting issue in `prism_pruner` 0.1.2 that can drop a unique conformer) and always retains the lowest-energy member of a duplicate pair. The moment-of-inertia stage is opt-in (`prism_moi: false` by default).
- `Criterion` column in the duplicate-removal log table.

### Changed
- **Supported Python: 3.12 – 3.14** (`requires-python = ">=3.12,<3.15"`). Python 3.11 is no longer supported (`prism_pruner` requires ≥ 3.12); 3.14 is newly supported. Upper bound kept below 3.15 because PyTorch and TBLite publish no 3.15 wheels yet. On 3.14 the `skala` extra resolves to `skala` 2026.8 (2026.9 declares `<3.14`). ReadTheDocs build moved to Ubuntu 24.04 / Python 3.12.

### Fixed
- Distance-matrix eigenvalues (clustering features and logged duplicate metric) are now computed with `numpy.linalg.eigvalsh` and sorted; `numpy.linalg.eig` returned them in arbitrary order, so feature columns did not correspond across conformers.
- Clustering uses `scipy.spatial.distance.cdist` instead of `scipy.spatial.distance_matrix` (deprecated in SciPy 1.18, removed in 1.20).
- `theory.md`: the distance matrix uses plain (not squared) distances; logged metric renamed `λ(D) RMSD [Å]`.
- Removed undocumented `-t` / `--threshold`, `--interest-vibro`, and `--interest-electro` CLI arguments from `_parser_arguments.py`.
- Synchronized missing protocol parameters (`fmax`, `maxiter`, `maxstep`, `validators`, `nroots`, `tda`, `read_population`) in `protocol.md`.
