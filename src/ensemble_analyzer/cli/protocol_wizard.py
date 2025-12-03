#!/usr/bin/env python3
"""
Interactive Protocol Wizard for Ensemble Analyser
Requires: Python >= 3.10, InquirerPy
"""

import json
import os

from InquirerPy import inquirer
from typing import Optional, Any
from importlib.resources import files


def load_grouped(path: str) -> list[str]:
    """Load grouped choices from file with category comments."""
    if not os.path.exists(path):
        return []
    
    grouped = []
    current_class = ""
    
    try:
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                
                # Category marker
                if line.startswith("#!"):
                    current_class = line[2:].strip()
                    continue
                
                # Skip comments and empty lines
                if not line or line.startswith(";") or line.startswith("#"):
                    continue
                
                # Process valid entries
                if current_class:
                    parts = line.split("|")
                    if len(parts) > 1:
                        formatted = f"{parts[0].strip():<30} # Defined for: {parts[1].strip():<40}"
                    else:
                        formatted = parts[0].strip()
                    grouped.append(f"{formatted:<40} # {current_class:<20}")
                else:
                    grouped.append(line)
    except Exception as e:
        print(f"⚠ Warning: Could not load {path}: {e}")
        return []
    
    return grouped


def load_functionals() -> list[str]:
    """Load available functionals from parameters file."""
    path = files("ensemble_analyzer").joinpath("parameters_file/functionals")
    return load_grouped(path)


def load_basis_sets() -> list[str]:
    """Load available basis sets from parameters file."""
    path = files("ensemble_analyzer").joinpath("parameters_file/basis_sets")
    return load_grouped(path)


def parse_choice(choice: str) -> str:
    """Extract clean value from formatted choice string."""
    return choice.split(" # ")[0].strip()


def get_int_input(message: str, default: str = "1") -> int:
    """Safe integer input with validation and retry."""
    while True:
        try:
            value = inquirer.text(message=message, default=default).execute()
            return int(value)
        except ValueError:
            print(f"⚠ Error: '{value}' is not a valid integer. Please retry.")


def get_float_input(message: str, default: str = "", allow_empty: bool = True) -> Optional[float]:
    """Safe float input with optional empty value."""
    value = inquirer.text(message=message, default=default).execute().strip()
    
    if allow_empty and not value:
        return None
    
    try:
        return float(value)
    except ValueError:
        print(f"⚠ Error: '{value}' is not a valid number.")
        if allow_empty:
            return None
        return float(default) if default else 1.0


def get_internal_coordinates() -> list[list[int]]:
    """Prompt for internal coordinates to monitor (bond/angle/dihedral)."""
    coordinates = []
    
    while inquirer.confirm(
        message="Add internal coordinate to monitor?",
        default=False
    ).execute():
        coord_input = inquirer.text(
            message="Coordinate (comma-separated atom indices, e.g., '1,2' or '1,2,3,4'):",
            validate=lambda x: len(x.strip()) > 0 and all(i.strip().isdigit() for i in x.split(','))
        ).execute()
        
        try:
            coord_list = [int(i.strip()) for i in coord_input.split(',')]
            
            # Validate: 2-4 atoms (bond, angle, dihedral)
            if 2 <= len(coord_list) <= 4:
                coordinates.append(coord_list)
                coord_type = {2: "Bond", 3: "Angle", 4: "Dihedral"}[len(coord_list)]
                print(f"✓ Added {coord_type}: {coord_list}")
            else:
                print("⚠ Error: Must specify 2-4 atoms (bond/angle/dihedral)")
        except ValueError:
            print("⚠ Error: Invalid format. Use comma-separated integers.")
    
    return coordinates


def get_constrains() -> list[list[int]]:
    """Prompt for geometric constraints."""
    constraints = []
    
    while inquirer.confirm(
        message="Add geometric constraint?",
        default=False
    ).execute():
        constraint_input = inquirer.text(
            message="Constraint (comma-separated atom indices, e.g., '1,2' or '1,2,3,4'):",
            validate=lambda x: len(x.strip()) > 0 and all(i.strip().isdigit() for i in x.split(','))
        ).execute()
        
        try:
            constraint_list = [int(i.strip()) for i in constraint_input.split(',')]
            
            if 2 <= len(constraint_list) <= 4:
                constraints.append(constraint_list)
                constraint_type = {2: "Bond", 3: "Angle", 4: "Dihedral"}[len(constraint_list)]
                print(f"✓ Added {constraint_type} constraint: {constraint_list}")
            else:
                print("⚠ Error: Must specify 2-4 atoms")
        except ValueError:
            print("⚠ Error: Invalid format. Use comma-separated integers.")
    
    return constraints


def clean_protocol_dict(step: dict[str, Any]) -> dict[str, Any]:
    """Remove None, empty values, and defaults from protocol dictionary."""
    cleaned = {}
    
    # Default values to skip
    skip_defaults = {
        'mult': 1,
        'charge': 0,
        'freq_fact': 1.0,
        'maxstep': 0.2,
        'calculator': 'orca',
        'graph': False,
        'no_prune': False,
        'cluster': False,
        'opt': False,
        'freq': False,
    }
    
    for key, value in step.items():
        # Skip None values
        if value is None:
            continue
        
        # Skip empty strings
        if isinstance(value, str) and not value:
            continue
        
        # Skip empty collections
        if isinstance(value, (list, dict)) and not value:
            continue
        
        # Skip default values
        if key in skip_defaults and value == skip_defaults[key]:
            continue
        
        cleaned[key] = value
    
    return cleaned


def print_step_summary(idx: int, step: dict[str, Any]) -> None:
    """Print formatted summary of protocol step."""
    print(f"\n{'='*60}")
    print(f"  STEP {idx} SUMMARY")
    print(f"{'='*60}")
    
    for key, value in step.items():
        if isinstance(value, dict):
            print(f"  {key}:")
            for k, v in value.items():
                print(f"    {k}: {v}")
        elif isinstance(value, list) and value and isinstance(value[0], list):
            print(f"  {key}:")
            for item in value:
                print(f"    {item}")
        else:
            print(f"  {key}: {value}")
    
    print(f"{'='*60}\n")


def protocol_step(step_num: int, level: str = "Basic") -> dict[str, Any]:
    """
    Create a single protocol step with level-dependent options.
    
    Args:
        step_num: Step number in the protocol
        level: Configuration level (Basic/Intermediate/Advanced)
    
    Returns:
        Dictionary containing all protocol parameters
    """
    step = {}
    
    print(f"\n{'#'*60}")
    print(f"  CONFIGURING STEP {step_num} - {level.upper()} MODE")
    print(f"{'#'*60}\n")
    
    # ==========================================
    # CALCULATOR SELECTION
    # ==========================================
    calculator = inquirer.select(
        message="Select calculator:",
        choices=["orca", "gaussian"],
        default="orca"
    ).execute()
    step["calculator"] = calculator
    
    # ==========================================
    # FUNCTIONAL & BASIS SET
    # ==========================================
    if calculator == "orca":
        # ORCA: Use fuzzy search from files
        functionals = load_functionals()
        if functionals:
            functional_raw = inquirer.fuzzy(
                message="Functional:",
                choices=functionals,
                validate=lambda x: len(x) > 0,
                instruction="Type to search or enter custom functional"
            ).execute()
            functional = parse_choice(functional_raw)
        else:
            functional = inquirer.text(
                message="Functional:",
                default="B3LYP"
            ).execute()
        
        # Dispersion correction
        if inquirer.confirm(
            message="Add dispersion correction?",
            default=False
        ).execute():
            dispersion = inquirer.select(
                message="Dispersion type:",
                choices=["D3", "D3BJ", "D4"],
                default="D4"
            ).execute()
            functional = f"{functional} {dispersion}"
        
        step["functional"] = functional
        
        # Basis set
        basis_sets = load_basis_sets()
        if basis_sets:
            basis_raw = inquirer.fuzzy(
                message="Basis set:",
                choices=basis_sets,
                validate=lambda x: len(x) > 0,
                instruction="Type to search or enter custom basis set"
            ).execute()
            step["basis"] = parse_choice(basis_raw)
        else:
            step["basis"] = inquirer.text(
                message="Basis set:",
                default="def2-SVP"
            ).execute()
    
    else:
        # Other calculators: Free text input
        step["functional"] = inquirer.text(
            message="Method/Functional:",
            default="B3LYP"
        ).execute()
        
        if inquirer.confirm(
            message="Add dispersion correction?",
            default=False
        ).execute():
            dispersion = inquirer.text(
                message="Input for dispersion - consider all keywords needed (e.g., D3, D3BJ, D4):",
                default="D3BJ"
            ).execute()
            step["functional"] = f"{step['functional']} {dispersion}"
        
        step["basis"] = inquirer.text(
            message="Basis set (basis-file is unsupported):",
            default=""
        ).execute()
    
    # ==========================================
    # SOLVENT
    # ==========================================
    if inquirer.confirm(
        message="Specify solvent?",
        default=False
    ).execute():
        solvent_dict = {}
        solvent_dict["solvent"] = inquirer.text(
            message="Solvent name:",
            default=""
        ).execute()
        solvent_dict["smd"] = inquirer.confirm(
            message="Use SMD solvation model?",
            default=False
        ).execute()
        step["solvent"] = solvent_dict
    
    # ==========================================
    # CALCULATION TYPE
    # ==========================================
    step["opt"] = inquirer.confirm(
        message="Perform geometry optimization?",
        default=False
    ).execute()
    
    step["freq"] = inquirer.confirm(
        message="Perform frequency calculation?",
        default=False
    ).execute()
    
    # ==========================================
    # BASIC PARAMETERS
    # ==========================================
    step["mult"] = get_int_input("Multiplicity:", default="1")
    step["charge"] = get_int_input("Charge:", default="0")
    
    read_orb = inquirer.text(
        message="Read orbitals from step (leave blank to skip):",
        default=""
    ).execute().strip()
    if read_orb:
        step["read_orbitals"] = read_orb
    
    add_input = inquirer.text(
        message="Additional calculator input (leave blank if none):",
        default=""
    ).execute().strip()
    if add_input:
        step["add_input"] = add_input.replace("\\n", "\n")
    
    # ==========================================
    # CLUSTERING & PRUNING
    # ==========================================
    if inquirer.confirm(
        message="Enable clustering?",
        default=False
    ).execute():
        cluster_num = get_int_input("Number of clusters:", default="5")
        step["cluster"] = cluster_num
    
    step["no_prune"] = inquirer.confirm(
        message="Disable pruning?",
        default=False
    ).execute()
    
    comment = inquirer.text(
        message="Comment (optional):",
        default=""
    ).execute().strip()
    if comment:
        step["comment"] = comment
    
    # ==========================================
    # INTERMEDIATE LEVEL
    # ==========================================
    if level in ["Intermediate", "Advanced"]:
        print("\n--- Intermediate Options ---")
        
        # Constraints
        constraints = get_constrains()
        if constraints:
            step["constrains"] = constraints
        
        # Monitor internals
        monitor = get_internal_coordinates()
        if monitor:
            step["monitor_internals"] = monitor
        
        # Skip optimization failures
        step["skip_opt_fail"] = inquirer.confirm(
            message="Skip failed optimizations?",
            default=False
        ).execute()
    
    # ==========================================
    # ADVANCED LEVEL
    # ==========================================
    if level == "Advanced":
        print("\n--- Advanced Options ---")
        
        thrG = get_float_input(
            "Energy threshold thrG in kcal/mol (leave blank for default):",
            allow_empty=True
        )
        if thrG is not None:
            step["thrG"] = thrG
        
        thrB = get_float_input(
            "Frequency threshold thrB in cm⁻¹ (leave blank for default):",
            allow_empty=True
        )
        if thrB is not None:
            step["thrB"] = thrB
        
        thrGMAX = get_float_input(
            "Maximum energy threshold thrGMAX in kcal/mol (leave blank for default):",
            allow_empty=True
        )
        if thrGMAX is not None:
            step["thrGMAX"] = thrGMAX
        
        step["freq_fact"] = get_float_input(
            "Frequency scaling factor (default 1.0):",
            default="1.0",
            allow_empty=False
        )
        
        read_pop = inquirer.text(
            message="Read population from step (leave blank to skip):",
            default=""
        ).execute().strip()
        if read_pop:
            step["read_population"] = read_pop

        skip_retention_rate = inquirer.confirm(
            message="Skip check on Minimum Retention Rate?",
            default=False
        ).execute()
        if skip_retention_rate: 
            step["skip_retention_rate"] = skip_retention_rate
    
    # Clean and return
    cleaned_step = clean_protocol_dict(step)
    print_step_summary(step_num, cleaned_step)
    
    return cleaned_step


def main():
    """Main entry point for protocol wizard."""
    print("="*60)
    print("  Ensemble Analyser Protocol Wizard")
    print("="*60)
    
    level = inquirer.select(
        message="Select configuration level:",
        choices=[
            "Basic - Essential parameters only",
            "Intermediate - Add constraints and monitoring",
            "Advanced - Full control over thresholds"
        ],
        default="Basic - Essential parameters only"
    ).execute().split(" - ")[0]
    
    protocol = {}
    step_num = 0
    
    while True:
        step = protocol_step(step_num, level=level)
        protocol[step_num] = step
        
        if not inquirer.confirm(
            message="Add another step to the protocol?",
            default=False
        ).execute():
            break
        
        step_num += 1
    
    # Save protocol
    filename = inquirer.text(
        message="Protocol filename:",
        default="protocol.json"
    ).execute()
    
    try:
        with open(filename, "w", encoding='utf-8') as f:
            json.dump(protocol, f, indent=4, ensure_ascii=False)
        print(f"\n✓ Protocol successfully saved to: {filename}")
        print(f"  Total steps: {len(protocol)}")
    except Exception as e:
        print(f"\n✗ Error saving protocol: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())