# Example 2: VCD and ECD Spectra Simulation

In this example, we simulate the Vibrational Circular Dichroism (VCD) and the Electronic Circular Dichroism (ECD) spectra of a chiral molecule and automatically fit it against an experimental reference file (`vcd_ref.dat`, `ecd_ref.dat`).

::::{tab-set}

:::{tab-item} Command
:sync: cmd

```bash
ensemble_analyzer -e ensemble.xyz -p protocol_vcd.json -cpu 44
```
:::

:::{tab-item} Protocol
:sync: proto

```json
{
    "0": {
        "functional": "B97-3c",
        "basis": "def2-mTZVP"
    },
    "1": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
        "opt": true,
        "freq": true
    },
    "2": {
        "functional": "wB97X-D4rev",
        "basis": "def2-QZVPP"
    },
    "3": {
        "functional": "B3LYP D4",
        "basis": "def2-QZVPPD",
        "freq": true
    }
}
```
:::

:::{tab-item} Ensemble
```{literalinclude} ../_static/ensemble_vcd.xyz
    :language: text
```
:::
:::{tab-item} Output
```{literalinclude} ../_static/output_vcd.out
    :language: text
```
:::

::::

At the end of the calculation, two files per spectra will be stored: one `.pickle` file and one `.png` file. If you want to modify the automatically generated spectra, it is possible to use the terminal interface of the `enan_graph_editor` module, to customize the spectra. 

<img src=../_static/vcd.png>

To refine the interested area and the multiplier, it is possible to change the values stored in the `setting.json` file and use the CLI command `enan_regraph` to recompute the various spectra, without the need of relaunching the calculation. 


The same workflow has been used to simulate the ECD spectra, but with the following ensemble and following protocol 

::::{tab-set}


:::{tab-item} Protocol
```json
{
    "0": {
	"functional": "r2SCAN-3c", 
	"basis": "def2-mTZVPP", 
	"opt": true, 
	"freq": true,
	"solvent":{
        "solvent": "acetonitrile", 
        "smd": true
	}},
    "1":{
        "functional": "wb97X-D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "acetonitrile",
            "smd": true
        },
        "add_input":"%tddft nroots 40 TDA FALSE end"
    },
    "2":{
        "functional": "CAM-B3LYP D4",
        "basis": "def2-TZVPP",
        "solvent": {
            "solvent": "acetonitrile",
            "smd": true
        },
        "add_input":"%tddft nroots 40 TDA false end",
        "read_orbitals": "1"
    }
}
```
:::

:::{tab-item} Ensemble
```{literalinclude} ../_static/ensemble_vcd.xyz
```
:::

::::

At the end of the calculations, the final ECD comparison figure is saved (both as a `.pickle` and a `.png`).

<img src=../_static/ecd.png>

