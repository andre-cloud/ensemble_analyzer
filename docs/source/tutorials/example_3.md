# Example 3: Singlet-Triplet Energy Gap (ΔE S→T)

This tutorial demonstrates how to calculate the **Singlet-Triplet energy gap** for an aminoborane derivative. The workflow involves manipulating the multiplicity settings across different steps to optimize Ground State (S$_0$), Excited Singlet (S$_1$), and Triplet (T$_1$) states.

::::{tab-set}

:::{tab-item} Command
```bash
ensemble_analyzer -e ensemble_aminoborane.xyz -p protocol_st_gap.json -cpu 16
```
:::

:::{tab-item} Ensemble
```{literalinclude} ../_static/ensemble_dest.xyz
    :language: text
```
:::

:::{tab-item} Protocol
```json
{
    "0": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "mult": 1,
        "charge": 0,
        "opt": true,
        "freq": true,
        "no_prune": true,
        "comment": "Opt Freq S0"
    },
    "1": {
        "functional": "CAM-B3LYP D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "mult": 1,
        "charge": 0,
        "no_prune": true,
        "comment": "SP S0"
    },
    "2": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "calculator": "orca",
        "mult": 1,
        "charge": 0,
        "opt": true,
        "freq": true,
        "read_orbitals": "2",
        "add_input": "\n! DeltaSCF UKS FreezeAndRelease SCFCheckGrad VerySlowConv\n %scf alphaconf 0,1 end",
        "no_prune": true,
        "comment": "Opt Freq S1"
    },
    "4": {
        "functional": "CAM-B3LYP D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "calculator": "orca",
        "mult": 1,
        "charge": 0,
        "opt": false,
        "freq": false,
        "no_prune": true,
        "comment": "Prepare SP S1"
    },
    "5": {
        "functional": "CAM-B3LYP D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "mult": 1,
        "charge": 0,
        "read_orbitals": "4",
        "add_input": "\n! DeltaSCF UKS FreezeAndRelease SCFCheckGrad\n %scf alphaconf 0,1 end %pal nprocs_group 4 end",
        "no_prune": true,
        "comment": "SP S1"
    },
    "6": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "mult": 3,
        "charge": 0,
        "opt": true,
        "freq": true,
        "no_prune": true,
        "comment": "Opt Freq T1"
    },
    "7": {
        "functional": "CAM-B3LYP D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "CHCl3",
            "smd": false
        },
        "calculator": "orca",
        "mult": 3,
        "charge": 0,
        "no_prune": true,
        "comment": "SP T1"
    }
}
```
:::

:::{tab-item} Output
```{literalinclude} ../_static/output_dest.out
    :language: text
```

::::

The workflow results in the energy values summarized in the diagram below:
<img src=../_static/jablonski_diagram.png> 

The final results can be schematically obtained using the CLI command `enan_get_energy`, which can also perform some basic algebric calculation between the average energies of the various protocols: 
```bash 
enan_get_energy --sub 4 1 \ # S0-S1 adiabatic gap
    --sub 7 1 \ # S0-T1 adiabatic gap
    --sub 5 7  # T1-S1 adiabatic gap
```

| Gap | Energy |
| :--: | :--: |
| S$_0$-S$_1$ gap |  ΔG = 64.4 kcal/mol| 
| S$_0$-T$_1$ gap |  ΔG = 62.1 kcal/mol| 
| T$_1$-S$_1$ gap |  ΔG = 0.10 eV |
