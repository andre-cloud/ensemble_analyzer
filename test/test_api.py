import pytest
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch

import ensemble_analyzer as ea
from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.protocol import Protocol
from ase.atoms import Atoms

def test_api_imports():
    assert hasattr(ea, "run")
    assert hasattr(ea, "EnsembleResult")
    assert hasattr(ea, "Conformer")
    assert hasattr(ea, "Protocol")
    assert hasattr(ea, "load_ensemble")
    assert hasattr(ea, "save_ensemble")
    assert hasattr(ea, "CalculationConfig")
    assert hasattr(ea, "__version__")

def test_conformer_ase_interop():
    atoms = Atoms(symbols="H2O", positions=[[0, 0, 0], [0, 1, 0], [0, 0, 1]])
    conf = Conformer.from_ase(atoms, number=1, raw=True)
    assert conf.number == 1
    assert conf.atoms == ("H", "H", "O")
    assert np.allclose(conf.geom, atoms.positions)

    atoms_out = conf.to_ase()
    assert list(atoms_out.get_chemical_symbols()) == ["H", "H", "O"]
    assert np.allclose(atoms_out.positions, atoms.positions)

@patch("ensemble_analyzer.api.read_ensemble")
def test_load_ensemble(mock_read):
    mock_read.return_value = ["mocked"]
    res = ea.load_ensemble("test.xyz")
    mock_read.assert_called_once_with("test.xyz", log=None, raw=False)
    assert res == ["mocked"]

@patch("ensemble_analyzer.api.save_snapshot")
def test_save_ensemble(mock_save):
    ea.save_ensemble("out.xyz", ["mocked"])
    mock_save.assert_called_once_with("out.xyz", ["mocked"], log=None)

@patch("ensemble_analyzer._managers.calculator_orchestration.CalculationOrchestrator")
@patch("ensemble_analyzer.api.create_logger")
def test_run_programmatic(mock_logger, mock_orch):
    # Mocking orchestrator properties
    instance = mock_orch.return_value
    
    c = Conformer(1, np.array([[0,0,0]]), ("C",), raw=True)
    instance.conformers = [c]
    
    res = ea.run(
        ensemble=[c],
        protocol={"0": {"functional": "B3LYP"}},
        quiet=True
    )
    
    assert isinstance(res, ea.EnsembleResult)
    assert len(res.conformers) == 1
    assert res.conformers[0] == c
    assert res.best == c
    instance.run.assert_called_once()
