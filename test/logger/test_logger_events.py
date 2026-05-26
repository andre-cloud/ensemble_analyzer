import pytest
import numpy as np
from unittest.mock import MagicMock
from datetime import timedelta
from ensemble_analyzer._logger.logger import Logger
from ensemble_analyzer._title import title


class TestLoggerEvents:

    @pytest.fixture
    def log(self):
        l = Logger("test_logger")
        l.info = MagicMock()
        l.debug = MagicMock()
        l.warning = MagicMock()
        l.error = MagicMock()
        l.critical = MagicMock()
        return l

    def test_title_screen(self, log):
        log.title_screen()
        log.info.assert_called_with(title)

    def test_application_input_received(self, log):
        config = {
            "conformers": 5, "len_protocols": 2, "temperature": 298.15,
            "cpu": 4, "restart": False, "protocols": []
        }
        log.application_input_received(config)
        assert log.info.call_count > 1

    def test_application_input_received_with_restart(self, log):
        config = {
            "conformers": 5, "len_protocols": 2, "temperature": 298.15,
            "cpu": 4, "restart": True, "protocols": []
        }
        log.application_input_received(config)
        assert any("RESTART" in str(c) for c in log.info.call_args_list)

    def test_application_correct_end(self, log):
        log.application_correct_end(timedelta(seconds=10), 5)
        log.info.assert_any_call("Total elapsed time: 0:00:10")
        log.info.assert_any_call("Final conformers: 5")

    def test_calculation_start(self, log):
        log.calculation_start(1, 1, 1)
        assert f"calc_1_1" in log._timers

    def test_calculation_failure(self, log):
        log.calculation_failure(1, "Some error occurred")
        args, _ = log.error.call_args
        assert "CONF 001" in args[0]
        assert "Some error occurred" in args[0]

    def test_missing_previous_thermo(self, log):
        log.missing_previous_thermo(1)
        args, _ = log.warning.call_args
        assert "conformer 1" in args[0]

    def test_missing_param(self, log):
        log.missing_param("temperature", "using default")
        args, _ = log.warning.call_args
        assert "temperature not found" in args[0]

    def test_pruning_start(self, log):
        log.pruning_start(1, 50)
        assert "pruning_1" in log._timers
        log.debug.assert_any_call("Starting pruning for protocol 1")

    def test_pruning_summary(self, log):
        log._timers["pruning_1"] = 0.0
        log.pruning_summary(1, 50, 30, 20)
        found = False
        for args, _ in log.info.call_args_list:
            if "Deactivated: 20" in args[0]:
                found = True
                break
        assert found

    def test_skip_pruning(self, log):
        log.skip_pruning(1)
        args, _ = log.warning.call_args
        assert "Pruning skipped" in args[0]

    def test_pca_analysis(self, log):
        log.pca_analysis(50, 5, True, "pca.pdf")
        log.info.assert_any_call("PCA Analysis:")
        log.info.assert_any_call("  Conformers: 50")

    def test_pca_analysis_no_clusters(self, log):
        log.pca_analysis(50, None, False, "pca.pdf")
        log.info.assert_any_call("  Include H: False")

    def test_spectra_generation(self, log):
        log.spectra_generation("IR", "ir.xy")
        args, _ = log.debug.call_args
        assert "IR spectrum" in args[0]

    def test_checkpoint_saved(self, log):
        log.checkpoint_saved(10)
        args, _ = log.debug.call_args
        assert "Checkpoint saved" in args[0]

    def test_checkpoint_loaded(self, log):
        log.checkpoint_loaded(10, 2)
        log.info.assert_any_call("Checkpoint loaded: 10 conformers")
        log.info.assert_any_call("Resuming from protocol 2")

    def test_spectra_start(self, log):
        log.spectra_start(1)
        assert "spectra_1" in log._timers

    def test_spectra_end(self, log):
        log._timers["spectra_1"] = 0.0
        log.spectra_end(1)
        found = False
        for args, _ in log.info.call_args_list:
            if "Sprectra convolution completed" in args[0]:
                found = True
                break
        assert found

    def test_spectra_skip(self, log):
        log.spectra_skip("IR")
        args, _ = log.warning.call_args
        assert "No calculation of IR graphs" in args[0]

    def test_converter_str(self, log):
        assert log.converter_str(3.14159) == "3.14"
        assert log.converter_str(5) == "5"

    def test_spectra_result(self, log):
        log.spectra_result("IR", {"Shift": 1.0, "FWHM": 10.0}, "Done")
        assert log.info.call_count >= 2

    def test_critical_error(self, log):
        log.critical_error("FILE_NOT_FOUND", "missing.out", path="/tmp")
        log.critical.assert_any_call("Error Type: FILE_NOT_FOUND")
        log.critical.assert_any_call("Message: missing.out")

    def test_table(self, log):
        log.table("Results", [["a", 1.0], ["b", 2.0]], ["Name", "Value"])
        assert log.info.call_count >= 2

    def test_timer_stop_unknown(self, log):
        assert log._stop_timer("nonexistent") == 0.0
