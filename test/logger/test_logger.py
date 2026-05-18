"""
Tests for the Logger module.
Verifies formatting, timer tracking, and context managers.
"""

import pytest
from unittest.mock import MagicMock
from datetime import timedelta
from ensemble_analyzer._logger.logger import Logger

class TestLogger:
    
    @pytest.fixture
    def log(self):
        """Return a Logger instance."""
        l = Logger("test_logger")
        l.info = MagicMock()
        l.debug = MagicMock()
        l.warning = MagicMock()
        l.error = MagicMock()
        l.critical = MagicMock()
        return l

    def test_timers(self, log):
        """Test internal timer logic."""
        log._start_timer("test")
        import time
        time.sleep(0.01)
        elapsed = log._stop_timer("test")
        assert elapsed > 0.0
        assert log._stop_timer("unknown") == 0.0

    def test_protocol_lifecycle(self, log):
        """Test logging of protocol start/end."""
        log.protocol_start(1, "OPT", "B3LYP", "6-31G", 10)
        assert "protocol_1" in log._timers
        log.info.assert_called()
        
        log.protocol_end(1, 10, 0)
        assert "protocol_1" not in log._timers

    def test_calculation_success_format(self, log):
        """Test formatting of calculation success message."""
        import numpy as np
        # Case 1: No imag freq (array present, but none < 0)
        log.calculation_success(1, 1, -100.0, -100.0, 5.0, np.array([10.0, 20.0]))
        call_args = log.info.call_args[0][0]
        assert "E = -100.00000000" in call_args
        # The logger prints "Imag. Freq 0" if frequencies are present but none are imaginary
        assert "Imag. Freq 0" in call_args 
        
        # Case 2: Imag freq
        log.calculation_success(1, 1, -100.0, -100.0, 5.0, np.array([-50.0, 10.0]))
        call_args = log.info.call_args[0][0]
        assert "Imag. Freq 1 (-50.00)" in call_args

    def test_track_operation_context(self, log):
        """Test track_operation context manager."""
        with log.track_operation("Complex Task", param=1):
            pass
        
        assert log.debug.call_count >= 2
        assert "Starting Complex Task (param=1)" in log.debug.call_args_list[0][0][0]

    def test_track_operation_failure(self, log):
        """Test exception handling in track_operation."""
        with pytest.raises(ValueError):
            with log.track_operation("Failing Task"):
                raise ValueError("Oops")
        
        # Check that error was logged with correct partial message
        assert log.error.called
        args, _ = log.error.call_args
        assert "Operation 'Failing Task' failed after" in args[0]