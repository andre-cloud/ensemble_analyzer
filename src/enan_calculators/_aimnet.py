import os
from pathlib import Path
from ._models import get_models_dir

try:
    import torch
    from aimnet.calculators import AIMNet2ASE
except ImportError:
    torch = None
    AIMNet2ASE = None


def create_aimnet_calc(charge, mult, method, solvent=None):
    if AIMNet2ASE is None:
        raise ImportError(
            "aimnet module missing. Install via: pip install aimnet[ase]@git+https://github.com/isayevlab/aimnetcentral.git"
        )

    model_path = get_models_dir("aimnet", create=False) / method

    if not model_path.exists():
        model_path = Path(method)
        if not model_path.exists():
            raise FileNotFoundError(f"AIMNet model not found: {model_path}")

    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", 1)))

    return AIMNet2ASE(
        str(model_path),
        charge=charge,
        mult=mult,
    )
