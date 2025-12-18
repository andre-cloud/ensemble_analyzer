# Quickstart

This tutorial guides you through a standard workflow: refining a conformational ensemble of a small molecule.

## Prerequisites

Ensure you have the following files in your directory:

1. `ensemble.xyz`: The initial conformers (from crest/rdkit).
2. `protocol.json`: The calculation settings.

## Step 1: Verify Setup

Check if your environment is ready:

```bash
enan_check
```

## Step 2: define the protocol

Create a `protocol.json` file. You can either use the protocol wizard: 
```bash
enan_protocol_wizard
```
Or use this standard template for: i) Single Point pruning, ii) Optimization + Frequency, and iii) re-evaluating the Electronic Energy:

```json
{
    "0": {
        "functional": "B97-3c",
        "solv": {"solvent": "chcl3", "smd": true}
    },
    "1": {
        "functional": "r2scan-3c",
        "opt": true,
        "freq": true,
        "solv": {"solvent": "chcl3", "smd": false}
    },
    "2": {
        "functional": "wB97X-D4rev", 
        "basis" : "def2-QZVPPD", 
        "solv": {"solvent": "chcl3", "smd": false}
    }
}
```

## Step 3: Run the program

Launch the main executable: 
```bash
ensemble_analyzer -e ensemble.xyz -p protocol.json -cpu 8
```

## Step 4 Check Result
After completion, check the log file `output.out`