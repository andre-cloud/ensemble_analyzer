import os
from pathlib import Path


def get_models_dir(calculator_name: str, create: bool = True) -> Path:
    env_path = os.getenv("ENAN_MODELS_DIR")
    if not env_path:
        env_path = os.getenv("ROT8_MODELS_DIR")
    base_dir = Path(env_path) if env_path else Path.home() / ".ensemble_analyzer" / "models"
    models_dir = base_dir / calculator_name
    print(f'{models_dir = }')
    if create:
        models_dir.mkdir(parents=True, exist_ok=True)
    return models_dir
