
import logging
import sys
from typing import Optional


class ColoredFormatter(logging.Formatter):
    """Log formatter with optional ANSI color support."""

    COLORS = {
        logging.DEBUG: "\033[90m",
        logging.INFO: "",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[1;31m",
    }
    RESET = "\033[0m"

    def __init__(self, fmt: Optional[str] = None, datefmt: Optional[str] = None, use_colors: Optional[bool] = None) -> None:
        """Configure the formatter with optional color auto-detection.

        Args:
            fmt: Log message format string.
            datefmt: Date/time format string.
            use_colors: None to auto-detect (stderr TTY), True to force, False to disable.
        """
        super().__init__(fmt, datefmt)
        
        if use_colors is None:
            self.use_colors = sys.stderr.isatty()
        else:
            self.use_colors = use_colors

    def format(self, record: logging.LogRecord) -> str:
        """Format a log record, optionally wrapping with ANSI color codes."""
        msg = super().format(record)
        
        if not self.use_colors:
            return msg
        
        color = self.COLORS.get(record.levelno, self.RESET)
        return f"{color}{msg}{self.RESET}"
    
    