# Changelog

All notable changes to this project will be documented in this file.


## [Unreleased]

### Added
- **Transition State (TS) Optimizer & Mode Validation**: Automated fragment-based localization of imaginary frequencies to validate transition states. Includes automatic displacement along spurious imaginary modes and subsequent re-optimization.
- **TD-DFT Support**: Added full support for time-dependent density functional theory (TD-DFT) calculations, including Tamm-Dancoff Approximation (TDA) control (`tda`) and configuration of the number of excited roots (`nroots`).
- **Electronic Spectral Generation**: Implemented Gaussian line-shape convolution for UV-Vis and Electronic Circular Dichroism (ECD) spectra generation.
- CLI Graph Editor (`enan_graph_editor`) and Spectra Regrapher (`enan_regraph`) tools.
- **ML Potentials Support**: Integration with AIMNet, UMA, FairChem, and MACE via the ASE calculator wrapper.
- Automatic Boltzmann weighting and conformational pruning strategies based on Euclidean Distance Matrix eigenvalues.
- **User templates in the Protocol Wizard**: protocols saved as `.json` in `~/.ensemble_analyzer/template/` (or `ENAN_TEMPLATES_DIR`) appear under *User Templates* next to the built-in *System Templates*. The folder layout is free (sub-folders and symlinked group folders are shown as sub-sections), and the new *Save as template* button stores the current protocol there.

### Fixed
- Removed undocumented `-t` / `--threshold`, `--interest-vibro`, and `--interest-electro` CLI arguments from `_parser_arguments.py`.
- Synchronized missing protocol parameters (`fmax`, `maxiter`, `maxstep`, `validators`, `nroots`, `tda`, `read_population`) in `protocol.md`.
