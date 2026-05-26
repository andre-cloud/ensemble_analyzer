"""
Global Pytest configuration and fixtures.
Provides shared mocks for Logger, Conformer, and Protocol to be used across the test suite.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock
from ensemble_analyzer.conformer.energy_data import EnergyRecord, EnergyStore
from ensemble_analyzer.conformer.conformer import Conformer

@pytest.fixture
def mock_logger():
    return MagicMock()

@pytest.fixture
def mock_conformer():
    """
    Returns a MagicMock simulating a Conformer object.
    Includes dunder methods for sorting and numeric returns for math operations.
    """
    conf = MagicMock()
    conf.active = True
    conf.number = 1
    conf.atoms = ["C", "O"]
    conf.geom = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    conf.write_xyz.return_value = "C 0.0 0.0 0.0\nO 1.2 0.0 0.0"
    
    # Mock energies return values (must be float for numpy/math)
    conf.get_energy.return_value = -100.0
    
    conf.energies = MagicMock()
    conf.energies.__contains__.return_value = False
    conf.energies.get.return_value = -100.0
    conf.energies.get_energy.return_value = -100.0 
    
    # Data structure for properties
    # IMPORTANT: Set E, G, H, zpve to floats to avoid "TypeError: < Mock" or numpy errors
    energy_record = MagicMock()
    energy_record.E = -100.0
    energy_record.G = -100.0
    energy_record.H = -100.0
    energy_record.zpve = 0.0
    energy_record.B = 10.0 
    energy_record.m = 1.0
    
    conf.energies.last.return_value = energy_record
    conf.energies.__getitem__.return_value = energy_record
    
    # Log methods
    conf.create_log.return_value = ["1", "-100.0"]

    # Enable sorting: Compare based on get_energy() return value
    conf.__lt__ = lambda self, other: self.get_energy() < other.get_energy()
    
    return conf

@pytest.fixture
def mock_protocol():
    proto = MagicMock()
    proto.number = 1
    proto.functional = "B3LYP"
    proto.basis = "6-31G*"
    proto.calculator = "orca" 
    proto.clustering = False
    proto.cluster = None
    proto.read_population = False
    proto.block_on_retention_rate = True
    proto.verbal_internals.return_value = []
    
    # Return tuple (calculator_object, label_string)
    proto.get_calculator.return_value = (MagicMock(), "opt")
    
    proto.thrG = 1.0
    proto.thrB = 1.0
    proto.thrGMAX = 999.0
    
    return proto


@pytest.fixture
def energy_record():
    return EnergyRecord(
        E=-100.0, G=-99.0, H=-98.5, zpve=0.05,
        B=1.5, m=2.0, Pop=50.0, Erel=0.0, time=5.0,
        Freq=np.array([100.0, 200.0, 300.0])
    )


@pytest.fixture
def energy_store():
    store = EnergyStore()
    store.add(1, EnergyRecord(E=-100.0, G=-99.0, H=-98.5, zpve=0.05, Pop=50.0))
    store.add(2, EnergyRecord(E=-200.0, G=-198.0, Pop=80.0))
    return store


@pytest.fixture
def conformer():
    conf = Conformer(
        number=1,
        geom=np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]),
        atoms=("C", "O"),
        raw=True
    )
    conf.energies.add(1, EnergyRecord(E=-100.0, G=-99.0, Pop=50.0))
    conf.energies.add(2, EnergyRecord(E=-200.0, G=-198.0, Pop=50.0))
    return conf