import pytest
import logging
from unittest.mock import patch, MagicMock
from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer._logger.formatter import ColoredFormatter


class TestColoredFormatter:

    def test_no_color(self):
        fmt = ColoredFormatter("%(message)s", use_colors=False)
        record = logging.LogRecord("test", logging.INFO, "", 0, "hello", (), None)
        assert fmt.format(record) == "hello"

    def test_with_color(self):
        fmt = ColoredFormatter("%(message)s", use_colors=True)
        record = logging.LogRecord("test", logging.WARNING, "", 0, "hello", (), None)
        result = fmt.format(record)
        assert result.startswith("\033[33m")
        assert result.endswith("\033[0m")

    def test_color_by_level(self):
        fmt = ColoredFormatter("%(message)s", use_colors=True)
        for level, code in [(logging.DEBUG, "\033[90m"), (logging.WARNING, "\033[33m"),
                            (logging.ERROR, "\033[31m"), (logging.CRITICAL, "\033[1;31m")]:
            record = logging.LogRecord("test", level, "", 0, "msg", (), None)
            result = fmt.format(record)
            assert result.startswith(code), f"Level {level} should use {code}"

    def test_auto_detect_color(self):
        with patch("ensemble_analyzer._logger.formatter.sys.stderr.isatty", return_value=True):
            fmt = ColoredFormatter(use_colors=None)
            assert fmt.use_colors is True

        with patch("ensemble_analyzer._logger.formatter.sys.stderr.isatty", return_value=False):
            fmt = ColoredFormatter(use_colors=None)
            assert fmt.use_colors is False


class TestCreateLogger:

    @patch("ensemble_analyzer._logger.create_log.logging.FileHandler")
    @patch("ensemble_analyzer._logger.create_log.Logger")
    def test_create_logger_default(self, MockLogger, MockFileHandler):
        mock_instance = MockLogger.return_value
        result = create_logger("test.log")
        assert result is mock_instance
        MockLogger.assert_called_with(name="enan")
        MockFileHandler.assert_called_with("test.log", mode="a")

    @patch("ensemble_analyzer._logger.create_log.logging.FileHandler")
    @patch("ensemble_analyzer._logger.create_log.Logger")
    def test_create_logger_debug(self, MockLogger, MockFileHandler):
        result = create_logger("test.log", debug=True)
        handler = MockFileHandler.return_value
        handler.setLevel.assert_called_with(logging.DEBUG)

    @patch("ensemble_analyzer._logger.create_log.logging.FileHandler")
    @patch("ensemble_analyzer._logger.create_log.Logger")
    def test_create_logger_custom_name(self, MockLogger, MockFileHandler):
        create_logger("test.log", logger_name="custom")
        MockLogger.assert_called_with(name="custom")

    @patch("ensemble_analyzer._logger.create_log.logging.FileHandler")
    @patch("ensemble_analyzer._logger.create_log.Logger")
    def test_create_logger_no_color(self, MockLogger, MockFileHandler):
        create_logger("test.log", disable_color=True)
        MockLogger.return_value.addHandler.assert_called()
