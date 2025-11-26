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
    
    # Create logger instance
    log = Logger(name=logger_name)
    log.basicConfig(
        filename=output_file,
        level=logging.DEBUG if DEBUG else logging.INFO,
        format=LOG_FORMAT,
        filemode="w",
    )
    
    # Disable noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    
    log.debug(f"Logger initialized | Debug: {debug} | Output: {output_file}")
    
    return log
