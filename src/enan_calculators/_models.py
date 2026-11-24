import os
from pathlib import Path


def get_models_dir(calculator_name: str) -> Path:
    env_path = os.getenv("ENAN_MODELS_DIR") or os.getenv("ROT8_MODELS_DIR")
    base_dir = Path(env_path) if env_path else Path.home() / ".ensemble_analyzer" / "models"
    return base_dir / calculator_name
