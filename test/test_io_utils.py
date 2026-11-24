"""
Tests for IO Utilities.
Verifies JSON serialisation.
"""

import io
import json

import numpy as np
import pytest

from ensemble_analyzer.io_utils import write_json, _serialise
from ensemble_analyzer.protocol.solvent import Solvent

class TestIOUtils:

    def test_serialise_int_keys(self):
        """_serialise converts integer dict keys to strings."""
        data = {0: {"E": -1961.03, "G": -1960.48}, 1: {"E": -1960.12}}
        result = _serialise(data)
        assert "0" in result and "1" in result
        assert all(isinstance(k, str) for k in result)

    def test_serialise_np_int_keys(self):
        """_serialise converts numpy integer dict keys to strings."""
        data = {np.int64(0): {"E": -1961.03}}
        result = _serialise(data)
        assert "0" in result
        assert all(isinstance(k, str) for k in result)

    def test_write_json_int_keys(self):
        """write_json produces valid JSON with integer keys."""
        data = {0: {"E": -1961.03}, 2: {"E": -1960.48}}
        buf = io.StringIO()
        write_json(data, buf)
        result = json.loads(buf.getvalue())
        assert result["0"]["E"] == -1961.03
        assert result["2"]["E"] == -1960.48

    def test_write_json_dataclass(self):
        """write_json serialises dataclass objects like Solvent."""
        data = {
            "0": {
                "functional": "r2scan-3c",
                "solvent": Solvent("CHCl3"),
            }
        }
        buf = io.StringIO()
        write_json(data, buf)
        result = json.loads(buf.getvalue())
        assert result["0"]["solvent"] == {"solvent": "CHCl3", "smd": False}
        assert result["0"]["functional"] == "r2scan-3c"