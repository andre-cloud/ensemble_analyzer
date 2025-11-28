
import logging


class ColoredFormatter(logging.Formatter):
    COLORS = {
        logging.DEBUG: "\033[90m",
        logging.INFO: "\033[37m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[1;31m",
    }
    RESET = "\033[0m"

    def __init__(self, fmt=None, datefmt=None, use_colors=None):
        """
        Args:
            use_colors: None=auto-detect, True=force, False=disable
        """
        super().__init__(fmt, datefmt)
        
        if use_colors is None:
            self.use_colors = self._should_use_colors()
        else:
            self.use_colors = use_colors

    def format(self, record):
        msg = super().format(record)
        
        if not self.use_colors:
            return msg
        
        color = self.COLORS.get(record.levelno, self.RESET)
        return f"{color}{msg}{self.RESET}"
    
    