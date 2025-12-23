"""
Tests for Configuration classes.
"""
import pytest
from unittest.mock import MagicMock
from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer._clustering.cluster_config import ClusteringConfig

def test_calculation_config_defaults():
    config = CalculationConfig()
    assert config.temperature == 298.15
    assert config.cpu > 0
    assert config.include_H is True

def test_protocol_calculator_retrieval():
    """Test the get_calculator helper on Protocol."""
    # Setup protocol with explicit parameters
    p = Protocol(number=1, opt=True, functional="B3LYP", basis="def2-SVP")
    
    # Verify property works based on opt=True
    assert p.calculation_level == "OPT"
    
    # FIX: Set an INVALID calculator name to force the lookup error logic
    # without trying to instantiate a real ASE calculator (which would crash if ORCA is missing)
    p.calculator = "fictitious_calc"
    
    conf = MagicMock()
    
    # Test that requesting a non-existent calculator raises ValueError
    # (Checking the logic in src/ensemble_analyzer/_protocol/protocol.py around line 220)
    with pytest.raises(ValueError) as exc: 
        p.get_calculator(1, conf)
    
    assert "not yet registered" in str(exc.value)

def test_clustering_config():
    c = ClusteringConfig(n_clusters=5)
    assert c.n_clusters == 5
    assert c.set_cluster_attribute is True