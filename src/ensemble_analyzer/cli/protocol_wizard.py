#!/opt/miniconda3/bin/python
 
from InquirerPy import inquirer
from InquirerPy.separator import Separator
import json, os

def load_grouped(path):
    if not os.path.exists(path):
        return []
    grouped = []
    current_class = ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#!"):
                current_class = line[2:].strip()
            if not line or line.startswith(";") or line.startswith("#"):
                continue
            elif not line.startswith("#"):
                # Aggiungi la categoria come prefisso, separata da " | "
                if current_class:
                    l = line.split("|")
                    if len(l) > 1:
                        line = f"{l[0].strip():<30} # Defined for: {l[1].strip():<40}"
                    else:
                        line = l[0].strip()
                    grouped.append(f"{line:<40} # {current_class:<20}")
                else:
                    grouped.append(line)
    return grouped

def load_functionals():
    path = os.path.join(os.path.dirname(__file__), "src", "parameters_file", "functionals")
    return load_grouped(path)

def load_basis_sets():
    path = os.path.join(os.path.dirname(__file__), "src", "parameters_file", "basis_sets")
    return load_grouped(path)    

def print_step_summary(idx, step):
    print(f"\n--- Step {idx} summary ---")
    for k, v in step.items():
        print(f"{k}: {v}")
    print("-------------------\n")

def protocol_step(step_num, level="Basic"):
    step = {}
    # BASIC
    functionals = load_functionals()
    if functionals:
        step["functional"] = inquirer.fuzzy(
            message="Functional:",
            choices=functionals,
            multiselect=False,
            validate=lambda x: len(x) > 0,
            instruction="Type to search or enter a new functional"
        ).execute().split(" # ")[0].strip()

    basis_sets = load_basis_sets()
    if basis_sets:
        step["basis"] = inquirer.fuzzy(
            message="Basis set:",
            choices=basis_sets,
            validate=lambda x: len(x) > 0,
            instruction="Type to search or enter a new basis set"
        ).execute().split(" # ")[0].strip()

    if inquirer.confirm(message="Do you want to specify a solvent?", default=False).execute():
        solv = {}
        solv["solvent"] = inquirer.text(message="Solvent name:", default="").execute()
        solv["smd"] = inquirer.confirm(message="Use SMD?", default=False).execute()
        step["solvent"] = solv
    else:
        step["solvent"] = {}
    step["opt"] = inquirer.confirm(message="Optimization using ASE optimizer?", default=False).execute()
    step["freq"] = inquirer.confirm(message="Frequency calculation?", default=False).execute()
    step["add_input"] = inquirer.text(message="Other ORCA input (leave blank if not needed):", default="").execute().replace("\\n", "\n")
    step["mult"] = int(inquirer.text(message="Multiplicity:", default="1").execute())
    step["charge"] = int(inquirer.text(message="Charge:", default="0").execute())
    step["read_orbitals"] = inquirer.text(message="Protocol number to read orbitals from (optional):", default="").execute()
    # cluster: False o integer
    use_cluster = inquirer.confirm(message="Clustering?", default=False).execute()
    if use_cluster:
        cluster_num = inquirer.text(message="How many clusters?", default="2").execute()
        try:
            step["cluster"] = int(cluster_num)
        except ValueError:
            step["cluster"] = 5
    else:
        step["cluster"] = False
    step["no_prune"] = inquirer.confirm(message="Disable pruning?", default=False).execute()
    step["comment"] = inquirer.text(message="Comment (optional):", default="").execute()

    # INTERMEDIATE
    if level in ["Intermediate", "Advanced"]:
        step["constrains"] = []
    else:
        step["constrains"] = []

    # ADVANCED
    if level == "Advanced":
        step["maxstep"] = float(inquirer.text(message="maxstep (default 0.2):", default="0.2").execute())
        thrG = inquirer.text(message="thrG (kcal/mol, leave blank for default):", default="").execute()
        step["thrG"] = float(thrG) if thrG else None
        thrB = inquirer.text(message="thrB (cm-1, leave blank for default):", default="").execute()
        step["thrB"] = float(thrB) if thrB else None
        thrGMAX = inquirer.text(message="thrGMAX (kcal/mol, leave blank for default):", default="").execute()
        step["thrGMAX"] = float(thrGMAX) if thrGMAX else None
        step["freq_fact"] = float(inquirer.text(message="Frequency scaling factor:", default="1.0").execute())
    else:
        step["thrG"] = None
        step["thrB"] = None
        step["thrGMAX"] = None
        step["freq_fact"] = 1.0
        step["maxstep"] = 0.2

    step["graph"] = False  # puoi aggiungere la domanda se vuoi
    step["calculator"] = 'orca'

    print_step_summary(step_num, step)
    return step

def main():
    protocol = {}
    step_num = 0
    print("Interactive wizard for Ensemble_Analyser protocol file creation")
    level = inquirer.select(
        message="Select configuration level:",
        choices=["Basic", "Intermediate", "Advanced"],
        default="Basic"
    ).execute()
    while True:
        step = protocol_step(step_num, level=level)
        protocol[step_num] = step
        another = inquirer.confirm(message="Add another step?", default=False).execute()
        if not another:
            break
        step_num += 1

    filename = inquirer.text(message="Protocol file name to save:", default="protocol.json").execute()
    with open(filename, "w") as f:
        json.dump(protocol, f, indent=4)
    print(f"Protocol saved in {filename}")

if __name__ == "__main__":
    main()