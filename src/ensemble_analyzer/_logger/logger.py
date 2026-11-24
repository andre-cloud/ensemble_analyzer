import logging
from tabulate import tabulate
import time
from typing import Generator, Optional, Any, Dict, List, Union
from contextlib import contextmanager
from datetime import timedelta



from ensemble_analyzer.constants import DEBUG, LOG_FORMAT



level_default = logging.DEBUG if DEBUG else logging.INFO



class Logger(logging.Logger): 
    """
    Custom Logger class with enhanced formatting and specialized event methods for EnAn.
    """

    ARROW = "→"
    TICK = "✓"
    FAIL = "✗"
    WARNING = "⊗"
    SPLIT = "|"


    def __init__(self, name: str, level: int = level_default) -> None:
        """Initialize the custom Logger.

        Args:
            name: Logger name.
            level: Logging level. Defaults to DEBUG or INFO based on settings.
        """
        super().__init__(name, level=level)

        self._timers: Dict[str, float] = {}
    # ===
    # Application Event
    # ===
    
    def application_input_received(self, config: Dict[str,Any]) -> None: 
        """
        Log the initial configuration summary.

        Args:
            config (Dict[str, Any]): Dictionary of run parameters (ensemble size, protocols, etc.).

        Returns:
            None
        """

        self._separator("Calculation Input")
        self.info(f"Ensemble: {config.get('conformers', 'N/A')} confromer(s)")
        self.info(f"Protocols: {config.get('len_protocols', 'N/A')}")
        self.info(f"Temperature: {config.get('temperature', 'N/A')} K")
        self.info(f"CPU cores: {config.get('cpu', 'N/A')}")
        if config.get('restart', False):
            self.info("Mode: RESTART")
        self.info('')
        self.info('Protocol Steps:')
        for protocol in config.get('protocols', []):
            self.info(f"{protocol.number:03d}. {self.ARROW} {str(protocol)} {self.SPLIT} {protocol.calculation_level} " +
                      (f"{self.SPLIT} {protocol.comment}" if protocol.comment else "" ) + 
                      f"\n {protocol.thr}")            
        self._separator()

    def application_correct_end(self, total_time: timedelta, total_conformers: int) -> None:
        """Log the successful end of the application with timing."""
        self._separator("Calculation COMPLETED")
        self.info(f"Total elapsed time: {total_time}")
        self.info(f"Final conformers: {total_conformers}")
        self._separator()


    # ===
    # Protocols Event
    # ===
    
    def protocol_start(self, number: int, level: str, functional:str, basis: str, active_conformers:int) -> None:
        """
        Log the start of a protocol step.

        Args:
            number (int): Protocol ID.
            level (str): Calculation level (e.g. "OPT+FREQ").
            functional (str): DFT functional.
            basis (str): Basis set.
            active_conformers (int): Number of conformers entering this step.

        Returns:
            None
        """
        
        self._separator(f"PROTOCOL {number} - {level}")
        self.info(f"Level: {functional}/{basis}")
        self.info(f"Active conformers: {active_conformers}")
        self._separator()
        self._start_timer(f"protocol_{number}")

    def protocol_end(self, number: int, active_conformers: int, deactivated: int) -> None:
        """Log the end of a protocol step with timing and counts."""
        elapsed = self._stop_timer(f"protocol_{number}")
        self.info("")
        self.info(f"Protocol {number} completed in {timedelta(seconds=elapsed)}")
        self.info(f"Active: {active_conformers} {self.SPLIT} Deactivated: {deactivated}")
        self._separator(f"END PROTOCOL {number}")
        self.info("")

    # === 
    # Calculation Events
    # ===

    def calculation_start(self, conformer_id: int, protocol_number: int, count: int) -> None:
        """Log the start of a single conformer calculation."""
        self.info(f"{count:03d}. {self.ARROW} CONF {conformer_id:03d} {self.SPLIT} Protocol {protocol_number}")
    
    def calculation_success(self, conformer_id: int, protocol_number: int, energy: float, gibbs:float, elapsed_time: float, frequencies: List[float]) -> None:
        """
        Log a successful calculation result.

        Args:
            conformer_id (int): Conformer ID.
            protocol_number (int): Protocol ID.
            energy (float): Electronic energy.
            gibbs (float): Gibbs free energy.
            elapsed_time (float): Duration in seconds.
            frequencies (List[float]): List of frequencies (for imaginary check).

        Returns:
            None
        """

        text = f"\t{self.TICK} CONF {conformer_id:03d} {self.SPLIT} E = {energy:.8f} Eh {self.SPLIT} Time: {elapsed_time:.1f}s"
        if frequencies is not None and len(frequencies) > 0: 
            n_im_freq = len(frequencies[frequencies<0])
            im_freq = f'({", ".join([f"{i:.2f}" for i in frequencies[frequencies<0]])})' if n_im_freq > 0 else ""
            text += f' {self.SPLIT} Imag. Freq {n_im_freq} {im_freq}'
        self.info(text)
    
    def missing_previous_thermo(self, conformer_id: int) -> None:
        """Warn that thermochemical data is unavailable for a conformer."""
        self.warning(f'{self.WARNING} No previous thermochemical data found for conformer {conformer_id}: setting G, H, S, ZPVE to NaN.')

    
    
    # ===
    # Pruning Events
    # ===

    def pruning_start(self, protocol_number: int, conformer_count: int) -> None:
        """Log the start of the pruning stage."""
        self.info("")
        self._separator("Pruning sequence", width=30, char="~")
        self.debug(f"Starting pruning for protocol {protocol_number}")
        self.debug(f"Conformers before pruning: {conformer_count}\n")
        self._start_timer(f"pruning_{protocol_number}")

    def pruning_summary(self, protocol_number: int, initial_count: int, final_count: int, deactivated_count: int) -> None:
        """
        Log statistics after the pruning stage.

        Args:
            protocol_number (int): Protocol ID.
            initial_count (int): Conformers before pruning.
            final_count (int): Conformers after pruning.
            deactivated_count (int): Number of removed conformers.

        Returns:
            None
        """

        elapsed = self._stop_timer(f"pruning_{protocol_number}")
        retention = (final_count / initial_count * 100) if initial_count > 0 else 0
        
        self.info(f"\n{self.TICK} Pruning completed in {elapsed:.4f}s")
        self.info(f"Initial: {initial_count} {self.ARROW} Final: {final_count} ({retention:.1f}% retained)"
        )
        self.info(f"Deactivated: {deactivated_count}")
        self.info("")

    def skip_pruning(self, protocol_number: int) -> None:
        """Log that pruning was skipped for a protocol."""
        self.warning(f'{self.WARNING} Pruning skipped for Protocol Step {protocol_number}.')

    # ===
    # Analysis Events
    # ===

    def pca_analysis(self, conformer_count: int, n_clusters: Optional[int], include_hydrogen: bool, output_file: str) -> None:
        """Log PCA analysis parameters."""
        self.info("")
        self.info("PCA Analysis:")
        self.info(f"  Conformers: {conformer_count}")
        if n_clusters:
            self.info(f"  Clusters: {n_clusters}")
        self.info(f"  Include H: {include_hydrogen}")
        self.info(f"  Output: {output_file}")
    
    # ===
    # Checkpoint Events
    # ===
    
    def checkpoint_saved(self, conformer_count: int) -> None:
        """Log that a checkpoint was saved."""
        self.debug(f"Checkpoint saved: {conformer_count} conformers ")
    
    # ===
    # Spectra Events
    # ===

    def spectra_start(self, protocol_number: int) -> None:
        """Log the start of spectra convolution."""
        self._separator("Spectra convolution", width=40, char="-")
        self.debug(f"Starting spectra convolution for protocol {protocol_number}")
        self._start_timer(f"spectra_{protocol_number}")

    def spectra_end(self, protocol_number: int) -> None:
        """Log the end of spectra convolution with timing."""
        elapsed = self._stop_timer(f"spectra_{protocol_number}")
        self.info(f"\n{self.TICK} Sprectra convolution completed in {elapsed:.4f}s")
        self.info("")

    def spectra_skip(self, graph_type: str) -> None:
        """Log that a graph type will be skipped."""
        self.warning(f'{self.WARNING} No calculation of {graph_type} graphs. Skipping')

    def converter_str(self, i: Union[int, float]) -> str:
        """Format a numeric value for display."""
        if isinstance(i, float):
            return f"{i:.2f}"
        return f"{i}"

    def spectra_result(self, graph_type: str, parameters: Dict, msg: str) -> None:
        """Log the result of a spectra convolution with its parameters.

        Args:
            graph_type: Type of spectrum (IR, VCD, UV, ECD).
            parameters: Dictionary of convolution parameters.
            msg: Informational message about the result.
        """
        res = [f"{k}: {self.converter_str(v)}" for k, v in parameters.items()]
        self._separator(f"{graph_type} Spectra convolution", char="-", width=35)
        self.info(msg)
        self.info("\t".join(res))


    # ===
    # Performance Tracking
    # ===

    def _start_timer(self, key: str) -> None:
        """Start a named timer."""
        self._timers[key] = time.perf_counter()
    
    def _stop_timer(self, key: str) -> float:
        """Stop a named timer and return the elapsed time.

        Args:
            key: Timer identifier.

        Returns:
            float: Elapsed time in seconds, or 0.0 if timer not found.
        """
        if key not in self._timers:
            return 0.0
        elapsed = time.perf_counter() - self._timers[key]
        del self._timers[key]
        return elapsed

    # ===
    # Formatting
    # ===

    def _separator(self, title: str = "", char: str = "=", width: int = 70):
        """Print a separator line."""
        if title:
            self.info(f"\n{char * width}")
            self.info(f"{title:^{width}}")
            self.info(f"{char * width}")
        else:
            self.info(f"{char * width}")
    
    def table(self, title: str, data: List[Any], headers: Union[List[str], str], **kwargs)-> None:
        """
        Print a formatted ASCII table.

        Args:
            title (str): Title of the table.
            data (List[Any]): List of rows (lists or tuples).
            headers (Union[List[str], str]): List of column headers.
            **kwargs: Additional formatting arguments for the separator.

        Returns:
            None
        """
        self._separator(title, **kwargs)
        self.info("\n" + tabulate(data, headers=headers, floatfmt=".10f", disable_numparse=True))


    @contextmanager
    def track_operation(self, operation_name: str, **context: Any) -> Generator[None, None, None]:
        """Context manager to track operation duration."""
        start = time.perf_counter()
        context_str = ", ".join(f"{k}={v}" for k, v in context.items())
        context_suffix = f" ({context_str})" if context_str else ""
        self.debug(f"Starting {operation_name}{context_suffix}")

        try:
            yield
        except Exception as e:
            elapsed = time.perf_counter() - start
            self.error(
                f"Operation '{operation_name}' failed after {elapsed:.2f}s{context_suffix}: {e}"
            )
            raise
        else:
            elapsed = time.perf_counter() - start
            self.debug(f"Completed {operation_name}{context_suffix} in {elapsed:.2f}s")


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

    def __init__(self, fmt: Optional[str] = None, datefmt: Optional[str] = None, use_colors: bool = False) -> None:
        super().__init__(fmt, datefmt)
        self.use_colors = use_colors

    def format(self, record: logging.LogRecord) -> str:
        """Format a log record, optionally wrapping with ANSI color codes."""
        msg = super().format(record)
        if not self.use_colors:
            return msg
        return f"{self.COLORS.get(record.levelno, self.RESET)}{msg}{self.RESET}"


def create_logger(
    output_file: str,
    debug: bool = False,
    logger_name: str = "enan",
    disable_color: bool = False,
) -> Logger:
    """
    Create the EnAn logger writing to *output_file*.

    Args:
        output_file: Output log filename
        debug: Enable debug level logging
        logger_name: Logger name (default: "enan")
    """
    log = Logger(name=logger_name)

    handler = logging.FileHandler(output_file, mode="a")
    handler.setLevel(logging.DEBUG if debug else logging.INFO)
    handler.setFormatter(ColoredFormatter(LOG_FORMAT, use_colors=disable_color))

    log.setLevel(logging.DEBUG if debug else logging.INFO)
    log.addHandler(handler)

    # Disable noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)

    log.debug(f"Logger initialized | Debug: {debug} | Output: {output_file}")
    return log
