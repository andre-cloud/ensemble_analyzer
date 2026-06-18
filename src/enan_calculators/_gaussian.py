from ase.calculators.gaussian import Gaussian


def create_gaussian_calc(charge, mult, method, basis, solvent=None, cpu=1, add_input="", label=None):
    solv = ""
    if solvent:
        solv = f" SCRF=(CPCM,Solvent={solvent})"

    route = f"# {method}/{basis}{solv}"

    if add_input.strip():
        route += " " + add_input.strip()

    return Gaussian(
        label=label or "gaussian",
        output_type='N',
        mem=f"{cpu * 2}GB",
        extra=route,
        charge=charge,
        mult=mult,
        nprocshared=cpu,
    )
