from .logger import Logger
import logging
import sys
from src.constants import DEBUG, LOG_FORMAT


def create_logger(
    output_file: str,
    debug: bool = False,
    also_console: bool = True,
    logger_name: str = "ensemble_logger"
) -> Logger:
    """
    Create enhanced logger - DROP-IN REPLACEMENT for old create_log().
    
    Args:
        output_file: Output log filename
        debug: Enable debug level logging
        also_console: Also output to console
        logger_name: Logger name (default: "ensemble_logger")
    
    Returns:
        Logger instance (subclass of logging.Logger)
    
    Usage:
        log = create_logger("output.log", debug=True)
        log.protocol_start(1, "OPT", "B3LYP", "def2-SVP", 50)
        log.info("Standard logging still works!")
    """
    global DEBUG
    DEBUG = debug
    
    # Register our custom logger class
    logging.setLoggerClass(Logger)
    
    # Create logger instance
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    
    # Remove existing handlers to avoid duplication
    logger.handlers.clear()
    
    # Create handlers
    handlers = [logging.FileHandler(output_file, mode='w', encoding='utf-8')]
    if also_console:
        handlers.append(logging.StreamHandler(sys.stdout))
    
    # Configure handlers
    formatter = logging.Formatter(LOG_FORMAT)
    for handler in handlers:
        handler.setLevel(logging.DEBUG if debug else logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    # Disable propagation to root logger
    logger.propagate = False
    
    # Disable noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)
    
    logger.debug(f"Logger initialized | Debug: {debug} | Output: {output_file}")
    
    return logger
