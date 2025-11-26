from .logger import Logger
import logging
import sys
from src.constants import DEBUG, LOG_FORMAT


def create_logger(
    output_file: str,
    debug: bool = False,
    logger_name: str = "enan"
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
    global DEBUG
    DEBUG = debug
    
    # Create logger instance
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    
    # Configure handlers
    formatter = logging.Formatter(LOG_FORMAT)
    
    # Disable noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    
    logger.debug(f"Logger initialized | Debug: {debug} | Output: {output_file}")
    
    return logger
