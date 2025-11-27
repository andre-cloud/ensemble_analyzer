
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

    def format(self, record, disable_color):
        if not disable_color:
            color = self.COLORS.get(record.levelno, self.RESET)
        else: 
            color = ""
            
        msg = super().format(record)
        return f"{color}{msg}{self.RESET}"