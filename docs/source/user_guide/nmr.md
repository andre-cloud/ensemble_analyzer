# NMR: shifts, couplings and spectra

An NMR step computes isotropic shieldings and/or spin-spin couplings J on every conformer
(ORCA, Gaussian or NWChem), averages them over the ensemble and writes tables and spectra.
Chemical shifts come from **reference molecules computed by EnAn at the same level**.

```json
{
  "0": {"calculator": "orca", "functional": "r2SCAN-3c", "opt": true, "freq": true},
  "1": {"calculator": "orca", "functional": "PBE0", "basis": "pcSseg-2",
        "solvent": {"solvent": "chloroform", "smd": false},
        "nmr": {
          "shielding": true,
          "coupling": true,
          "nuclei": ["1H", "13C", "19F"],
          "references": [
            {"name": "TMS", "nuclei": ["1H", "13C"]},
            {"name": "CFCl3", "nuclei": ["19F"]}
          ],
          "decouple": {"13C": ["1H"]}
        }}
}
```

An NMR step is a single point: `opt`, `freq`, `ts`, `nroots` and `deltascf` are not allowed on
it, and the NMR keywords must not be repeated in `add_input`. ML potentials and TBLite cannot
compute NMR properties.

## Options

| Key | Default | Meaning |
|---|---|---|
| `shielding` | `false` | isotropic shieldings → shifts (needs `references`) |
| `coupling` | `false` | spin-spin couplings J (no reference needed) |
| `nuclei` | required | isotopes, e.g. `["1H", "13C", "15N", "2H"]`; one isotope per element |
| `atoms` | all | compute only these atoms (0-based, like every index in the protocol) |
| `references` | — | reference molecules (below); one reference per nucleus |
| `coupling_from` | auto | for a shielding-only step: NMR step whose J are used for the spectrum (auto = closest earlier step with `coupling`; `false` = singlets) |
| `equivalent` | `[]` | extra groups of equivalent atoms (0-based) |
| `decouple` | `{}` | `{"observed": ["partners"]}`, e.g. `{"13C": ["1H"]}` |
| `coupling_terms` | `"all"` | `"fc"`: Fermi-contact term only (cheaper) |
| `frequency` | `400` | spectrometer 1H frequency (MHz); X nuclei use `frequency·|γX/γH|` |
| `fwhm` | `1.0` | line width in Hz, a number or `{"1H": 0.5, "13C": 2}` |
| `opt_from` | auto | step whose level optimises the references |

## References

Each reference is `{"name": ..., "nuclei": [...]}` plus one geometry source, in this priority:

1. `"xyz": "file.xyz"` (path relative to the working directory, `charge`/`mult` optional);
2. a **bundled** reference, by name (case-insensitive, aliases accepted): `TMS` (1H, 13C, 29Si),
   `CFCl3` (19F), `H3PO4` (31P), `CH3NO2` and `NH3` (15N), `BF3OEt2` (11B), `Me2Se` (77Se),
   `Me4Sn` (119Sn);
3. `"smiles": "..."` through RDKit (`conda install -c conda-forge rdkit` or
   `pip install "ensemble-analyzer[nmr]"`).

Optional: `atoms` (0-based atoms of the reference used for σ_ref, default every atom of the
element) and `delta` (δ of the reference on the scale, for secondary references; default 0).

The geometry is only a starting point. The reference is **optimised at the level of the last
optimisation step before the NMR step** (TS, ΔSCF and TD-DFT optimisations excluded; ML and
xTB optimisations included). If there is none, it is optimised at the level of the NMR step
and a warning is written. Its shielding is then computed at the level of the NMR step:

δ = (σ_ref − σ) / (1 − σ_ref·10⁻⁶) + δ_ref

A second NMR step reuses the optimised reference and only computes the new shielding; a new
optimisation level re-optimises it starting from the last geometry. Everything is kept in
`nmr_references/<name>/` (`reference.json` is the cache), so a restart computes nothing twice.

## Ensemble averaging and equivalent nuclei

Shieldings and couplings are Boltzmann-averaged atom by atom (fast exchange) with the
populations of the step, or of `read_population`. Equivalent nuclei are then averaged: the
point group of the heavy atoms on the most populated conformer (filtered by the bond
topology), freely rotating XH3 groups (CH3, NH3⁺), and the `equivalent` groups (use them for
CF3 groups or for exchange that the conformers do not sample).

## Spectra

With `shielding`, EnAn writes one spectrum per nucleus, `NMR_<isotope>_p<N>_comp.xy`, and at
the end `NMR_<isotope>_comparison.png` with every NMR step:

- homonuclear couplings: exact spin Hamiltonian (second order, roof effect), large networks
  split into clusters;
- heteronuclear couplings **present in the calculation**: weak coupling weighted by the natural
  abundance of the partner (H–F and H–P splittings, 119Sn satellites; partners below 5 %
  abundance such as 13C are ignored, and so are 13C–13C couplings);
- `decouple` removes partners (e.g. 13C{1H});
- without couplings (no `coupling` and no earlier coupling step) the lines are singlets.

An experimental spectrum `nmr_<isotope>_ref.dat` (two columns: ppm, intensity, e.g.
`nmr_1h_ref.dat`) in the working directory is overlaid; the computed spectrum is shifted to
match it and the shift and similarity are reported.

## Output

- `nmr_p<N>.out`: references (σ_ref, optimisation level), atoms (σ, δ, group), signals
  (averaged over equivalent nuclei), couplings (|J| ≥ 0.1 Hz) and conformer weights;
- `nmr_p<N>.json`: the same data, for further processing;
- spectra and plots as above. `enan_regraph N` recomputes an NMR step from the checkpoint and
  the reference cache.

## Engines

| | ORCA | Gaussian | NWChem |
|---|---|---|---|
| keywords | `! NMR` + `%eprnmr` after the coordinates, with `ist` per element | `NMR=(GIAO,SpinSpin[,ReadAtoms])` | `property` block, `task dft property` |
| couplings | isotopes as requested | from the reduced K | needs `odft`, added automatically |
| atom subset | `Nuclei = i, j {…}` (1-based) | `ReadAtoms` (1-based) | atom lists (1-based) |

J is always reported for the isotopes in `nuclei`, whatever isotope the engine used.
