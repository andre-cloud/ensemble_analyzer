# Example 4: Transition State Ensemble Refinement

Ensemble Analyzer can also handle **Transition State (TS)** refinement. This example reproduces the refinement of a TS ensemble for a model Diels-Alder reaction (673 conformers generated via CREST).

The key features used here are: i) **clustering** to reduce the ensemble size, ii) **constrained optimization** to preserve the reaction center before the final TS optimization, and iii) the actual saddle point optimization via triggering the OptTS keyword. 

::::{tab-set}

:::{tab-item} Command
```bash
ensemble_analyzer -e ts_ensemble.xyz -p protocol_ts.json -cpu 44
```
:::

:::{tab-item} Ensemble
```{literalinclude} ../_static/ensemble_ts.xyz
    :language: text
```
:::

:::{tab-item} Protocol
```json
{
    "0": {
        "functional": "xTB",
        "cluster": 30,
	    "comment": "Clustering"
   },
    "1": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
        "opt": true,
        "comment": "Constrained Optimization",
        "constrains": [0,1,2,3]
    },
    "2": {
        "functional": "r2SCAN-3c",
        "basis": "def2-mTZVPP",
	    "freq": true,
        "add_input": "\n! OptTS \n%geom calc_hess true end\n",
        "read_orbitals": "1",
        "comment": "Berny"
    },
    "3": {
        "functional": "wB97X-3c",
        "comment": "SP"
    }
}
```
*Note: The atom indices in the constraints block must be adjusted to match the reacting atoms in your specific xyz file. Numbering starts at 0.*
:::

:::{tab-item} Output
```{literalinclude} ../_static/output_ts.out
    :language: text
```
:::
::::

The low-lying conformer after the refinement of the ensemble has an energy 4.84 kcal/mol **lower** than the original TS obtained from the Nudged Elastic Band. 

<img src=../_static/ts.png>
