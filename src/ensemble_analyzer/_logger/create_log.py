from .logger import Logger
from .formatter import ColoredFormatter
import logging
import sys
from ensemble_analyzer.constants import DEBUG, LOG_FORMAT


def create_logger(
    output_file: str,
    debug: bool = False,
    logger_name: str = "enan", 
    disable_color: bool = False,
) -> Logger:
    """
    Create enhanced logger
    
    Args:
        output_file: Output log filename
        debug: Enable debug level logging
        logger_name: Logger name (default: "enan")
    
    Returns:
        Logger instance (subclass of logging.Logger)
    """
    
    # Create logger instance
    log = Logger(name=logger_name)

    # File handler
    handler = logging.FileHandler(output_file, mode="w")
    handler.setLevel(logging.DEBUG if debug else logging.INFO)
    formatter = ColoredFormatter(LOG_FORMAT, use_colors=disable_color)
    handler.setFormatter(formatter)

    # Attach
    log.setLevel(logging.DEBUG if debug else logging.INFO)
    log.addHandler(handler)

    # Disable noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)

    log.debug(f"Logger initialized | Debug: {debug} | Output: {output_file}")
    return log
