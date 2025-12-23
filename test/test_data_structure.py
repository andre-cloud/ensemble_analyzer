"""
Tests for Data Structures (EnergyStore, SpectralStore).
Covers storage, retrieval, and serialization logic.
"""

import pytest
import numpy as np
from datetime import datetime
from ensemble_analyzer._conformer.energy_data import EnergyStore, EnergyRecord
from ensemble_analyzer._conformer.spectral_data import SpectralStore, SpectralRecord

class TestDataStructures:

    def test_energy_record_defaults(self):
        """Test EnergyRecord initialization and defaults."""
        now = datetime.now()
        rec = EnergyRecord(E=-100.0, time=now)

        assert rec.E == -100.0
        assert rec.G is np.nan # Default
        assert not rec.Freq

    def test_energy_store_operations(self):
        """Test adding, retrieving, and iterating EnergyStore."""
        store = EnergyStore()
        now = datetime.now()
        
        # Test Add
        rec1 = EnergyRecord(E=-100.0, time=now)
        store.add(1, rec1)
        
        assert store.__contains__(1)
        assert not store.__contains__(2)
        assert len(store.data) == 1
        
        # Test Get
        fetched = store[1]
        assert fetched.E == -100.0
        
        # Test Last
        rec2 = EnergyRecord(E=-200.0, time=now)
        store.add(2, rec2)
        assert store.last().E == -200.0
        
        # Test Get Energy Helper
        # Should return G if present, else E
        assert store.get_energy() == -200.0 # G is nan
        
        store.data[2].G = -100.1
        assert store.get_energy() == -100.1 # G is set

    def test_energy_store_log_info(self):
        """Test data extraction for logging."""
        store = EnergyStore()
        rec = EnergyRecord(
            E=-100.0, G=-100.1, B=10.0, m=1.0, 
            time=1.0,
            Erel=0.0, Pop=50.0, G_E=-0.1
        )
        store.add(1, rec)
        
        # Expected tuple: E, G-E, G, B, Erel, Pop, time
        info = store.log_info(1)
        assert info[0] == -100.0 # E
        assert np.isclose(info[1], -0.1) # G-E
        assert info[5] == '50.00' # Pop as a string

    def test_spectral_store_operations(self):
        """Test SpectralStore adding and checking."""
        store = SpectralStore()
        rec = SpectralRecord(X=np.array([1,2]), Y=np.array([0.1, 0.2]))
        
        store.add(1, "IR", rec)
        
        assert store.__has_graph_type__(1, "IR")
        assert not store.__has_graph_type__(1, "UV")
        assert not store.__has_graph_type__(2, "IR")
        
        # Test retrieval
        fetched = store.__getitem__(protocol_number=1, graph_type="IR")
        assert np.allclose(fetched.X, [1, 2])
        
        # Test magic methods
        assert store.__getitem__(protocol_number=1, graph_type="IR") == fetched

    def test_energy_record_property_defaults(self):
        """Test fallback defaults for missing properties."""
        rec = EnergyRecord(E=-100.0)
        # Assuming defaults like m=1, B=0 if not specified in __init__
        # Adjust based on actual dataclass defaults
        assert rec.m == None