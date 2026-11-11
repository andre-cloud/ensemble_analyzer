# ΔSCF excited states (ORCA, NWChem)

A protocol step with the `deltascf` key optimises an excited state as a single
non-aufbau determinant (orbital-optimised DFT). EnAn takes care of the parts
that make ΔSCF fragile: the starting orbitals, the convergence strategy, the
detection of a collapse to the ground state and the spin purification.

```json
{
  "1": {"functional": "CAM-B3LYP", "basis": "def2-TZVP", "dispersion": "D4",
        "nroots": 5, "tda": true, "comment": "TD-DFT reference"},
  "2": {"functional": "CAM-B3LYP", "basis": "def2-TZVP", "dispersion": "D4",
        "read_orbitals": "1", "deltascf": {"excitation": "auto", "root": 1},
        "comment": "ΔSCF S1"}
}
```

`opt` and `freq` work as in any other step (excited-state optimisation and frequencies).

## Options

| Key | Default | Meaning |
|-----|---------|---------|
| `excitation` | `"auto"` | `"auto"`: dominant transition of TD-DFT root `root`; or `"HOMO->LUMO"`, `"HOMO-n->LUMO+m"` |
| `root` | `1` | TD-DFT root used by `auto`, relative to the step multiplicity (see below) |
| `spin` | `"alpha"` | spin of the excited electron for explicit excitations (`auto` reads it from the TD-DFT output) |
| `spin_purify` | `true` | spin purification, method chosen automatically (see below) |
| `doublet_branch` | `"lower"` | doublet kept after the recoupling: `"lower"` or `"upper"` |
| `cation_guess` | `true` | last attempts start from the orbitals of the cation |
| `min_gap_ev` | `0.1` | minimum ΔSCF − reference gap [eV] (collapse check) |
| `max_attempts` | `null` | limit the convergence ladder (null = all attempts) |
| `alphaconf`, `betaconf` | `null` | ORCA only: raw occupation lists, override `excitation` |

Orbital indices are never written by hand: excitations are relative labels,
translated at run time for each conformer (ORCA `alphaconf`/`betaconf`,
NWChem `swap` with the occupied orbitals read from the reference output, so
ECPs are handled).

## Reference (#0)

Every ΔSCF step starts from a **ground-state + TD-DFT reference at the same level
of theory and geometry**, which provides the ground-state orbitals (initial
guess), E(GS) (collapse check), the composition of the TD-DFT roots (`auto`) and
their energies (state-identity check).

* `read_orbitals` pointing to an earlier TD-DFT step (`nroots`) at the same level
  (functional, basis, dispersion, solvent, charge, multiplicity): reused.
* No `read_orbitals`, or a reference at another level: the reference is computed
  for each conformer (TD-DFT, TDA, `iroot 0`).
* `read_orbitals` pointing to a non-TD-DFT step or to another ΔSCF step is an error:
  re-reading ΔSCF orbitals applies the excitation twice and collapses to S0
  (observed: `<S²>` = 0).

**Multiplicity.** ΔSCF keeps M_s: the reference is the lowest state of the step
multiplicity and roots are numbered from it.

| `mult` | Reference | `root: 1` |
|--------|-----------|-----------|
| 1 | S0 | S1 |
| 3 | T1 | T2 (T1 itself is a plain `mult: 3` step) |
| 2 | D0 | D1 |

`auto` warns when the dominant transition weighs less than 0.6: the state is
multiconfigurational and a single determinant describes it poorly.

## Convergence ladder

Attempts are tried in order until one is accepted. Orbitals always come from the
reference (or the cation), never from a failed attempt.

| # | Guess | ORCA ≥ 6.1 | NWChem |
|---|-------|------------|--------|
| 1 | reference | `DeltaSCF GMF`, `SOSCFUpdateSPOEst false` | `max_ovl` (MOM) |
| 2 | reference | GMF, `SOSCFMaxStep 0.1` | `imom` |
| 3 | reference | `FreezeAndRelease SlowConv` | MOM, `damp 70 ncydp 50` |
| 4 | reference | `FreezeAndRelease VerySlowConv`, PMOM, `SOSCFMaxStep 0.05` | MOM from the cation |
| 5 | cation | GMF | IMOM from the cation, damped |
| 6 | cation | FreezeAndRelease, very slow | — |

ORCA 6.0 (no GMF/FreezeAndRelease) uses MOM/PMOM with `SlowConv`/`VerySlowConv`;
ORCA < 6 is rejected. GMF always runs with `SOSCFUpdateSPOEst false`: with the
automatic saddle-point order an optimisation crashed (order reset to 0,
`NVecs<=0. This is a bug`). In NWChem `nodiis` made the SCF diverge and is not used.

**Acceptance.** An attempt is accepted when the run converged (ORCA exits with
code 0 also on errors: the output text is checked), the gap to the reference is
above `min_gap_ev` and, for open-shell singlets, `<S²>` ≥ 0.2. A gap closer to a
lower TD-DFT root than to the target one is reported as a possible collapse onto a
lower excited state (warning: ΔSCF and TD-DFT differ by up to ~0.4 eV).

## Spin purification

| Target determinant | Method | Extra calculations |
|--------------------|--------|--------------------|
| all unpaired electrons parallel (e.g. T2 from T1, D1 from D0 through the SOMO) | none, already pure | 0 |
| open-shell singlet (two open shells) | Yamaguchi AP: E_S = E_mix + (E_mix − E_T) `<S²>`_mix / (`<S²>`_T − `<S²>`_mix) | 1 triplet SP |
| doublet with three open shells | recoupling of the quartet and the three mixed determinants (Hait et al., *J. Chem. Phys.* **153**, 134108, 2020) | 3 SCF |
| four or more open shells | not supported: warning | 0 |

The purified energy is stored as E of the step. With `opt`, the geometry is the one
of the mixed state (no purified gradient) and the purification is done at the final
geometry; with `freq`, G = E_S + (G − E) of the mixed state.

## Validation (formaldehyde, PBE)

| | ORCA 6.1.1 | NWChem 7.3.1 |
|-|------------|--------------|
| S1 gap (def2-SVP) | 3.51 eV, `<S²>` 1.006 | 3.51 eV, `<S²>` 1.006 |
| E_S − E_mix (Yamaguchi) | +0.166 eV | +0.166 eV |
| T2 from T1 | 2.88 eV, `<S²>` 2.005 | — |
| Cation D(π→π*) recoupled doublets | 7.35 / 9.33 eV | 7.35 / 9.32 eV |
| `alphaconf` ↔ `swap` | `0,1` / `0,1,1` / `0,0,1` | `8 9` / `7 9` / `8 10` (same energies) |

Failures were forced (`MaxIter 3`) to check that the ladder moves on.
Water and formaldehyde converge at the first attempt: the ladder matters for harder
cases (strong charge transfer).

Writing `DeltaSCF` by hand in `add_input` (without the `deltascf` key) is still
possible; mixing it with the `deltascf` key is rejected.
