import logging
import sys
from tabulate import tabulate
import time
from pathlib import Path
from typing import Optional, Any, Dict, List, Union
from contextlib import contextmanager
from datetime import timedelta



from ensemble_analyzer.constants import DEBUG
from ensemble_analyzer.title import title



level_default = logging.DEBUG if DEBUG else logging.INFO



class Logger(logging.Logger): 

    ARROW = "→"
    TICK = "✓"
    FAIL = "✗"
    WARNING = "⊗"
    SPLIT = "|"


    def __init__(self, name:str, level=level_default):
        super().__init__(name, level=level)

        self._timers: Dict[str, float] = {}
    # ===
    # Application Event
    # ===
    
    def title_screen(self):
        self.info(title)

    def application_input_received(self, config: Dict[str,Any]): 
        self._separator("Calculation Input")
        self.info(f"Ensemble: {config.get('conformers', 'N/A')} confromer(s)")
        self.info(f'Protocols: {config.get('len_protocols', 'N/A')}')
        self.info(f"Temperature: {config.get('temperature', 'N/A')} K")
        self.info(f"CPU cores: {config.get('cpu', 'N/A')}")
        if config.get('restart', False):
            self.info("Mode: RESTART")
        self.info('')
        self.info('Protocol Steps:')
        for protocol in config.get('protocols', []):
            self.info(f"{int(protocol.number):03d}. {self.ARROW} {str(protocol)} {self.SPLIT} {protocol.calculation_level} " +
                      (f"{self.SPLIT} {protocol.comment}" if protocol.comment else "" ) + 
                      f"\n {protocol.thr}")            
        self._separator()

    def application_correct_end(self, total_time: timedelta, total_conformers: int):
        self._separator("Calculation COMPLETED")
        self.info(f"Total elapsed time: {total_time}")
        self.info(f"Final conformers: {total_conformers}")
        self._separator()


    # ===
    # Protocols Event
    # ===
    
    def protocol_start(self, number: int, level: str, functional:str, basis: str, active_conformers:int):
        self._separator(f"PROTOCOL {number} - {level}")
        self.info(f"Level: {functional}/{basis}")
        self.info(f"Active conformers: {active_conformers}")
        self._separator()
        self._start_timer(f"protocol_{number}")

    def protocol_end(self, number: int, active_conformers: int, deactivated: int):
        elapsed = self._stop_timer(f"protocol_{number}")
        self.info("")
        self.info(f"Protocol {number} completed in {timedelta(seconds=elapsed)}")
        self.info(f"Active: {active_conformers} {self.SPLIT} Deactivated: {deactivated}")
        self._separator(f"END PROTOCOL {number}")
        self.info("")

    # === 
    # Calculation Events
    # ===

    def calculation_start(self, conformer_id: int, protocol_number: int, count: int):
        self.info(f"{count:03d}. {self.ARROW} CONF {conformer_id:03d} {self.SPLIT} Protocol {protocol_number}")
        self._start_timer(f"calc_{conformer_id}_{protocol_number}")
    
    def calculation_success(self, conformer_id: int, protocol_number: int, energy: float, gibbs:float, elapsed_time: float, frequencies: List[float]):
        self._stop_timer(f"calc_{conformer_id}_{protocol_number}")
        text = f"\t{self.TICK} E = {energy:.8f} Eh {self.SPLIT} Time: {elapsed_time:.1f}s"
        if frequencies.size > 0: 
            n_im_freq = len(frequencies[frequencies<0])
            im_freq = f'({", ".join([f"{i:.2f}" for i in frequencies[frequencies<0]])})' if n_im_freq > 0 else ""
            text += f' {self.SPLIT} Imag. Freq {n_im_freq} {im_freq}'
        self.info(text)
    
    def calculation_failure(self, conformer_id: int, error: str):
        self.error(f"    {self.FAIL} CONF {conformer_id:03d} [FAILED] {self.SPLIT} Error: \n{error[:60]}")

    def missing_previous_thermo(self, conformer_id:int):
        self.warning(f'{self.WARNING} No previous thermochemical data found for conformer {conformer_id}: setting G, H, S, ZPVE to NaN.')

    def missing_param(self, param:str, action:str):
        self.warning(f'{self.WARNING} {param} not found. {action}')
    
    
    # ===
    # Pruning Events
    # ===

    def pruning_start(self, protocol_number: int, conformer_count: int):
        self.info("")
        self._separator("Pruning sequence", width=30, char="~")
        self.debug(f"Starting pruning for protocol {protocol_number}")
        self.debug(f"Conformers before pruning: {conformer_count}\n")
        self._start_timer(f"pruning_{protocol_number}")

    def pruning_summary(self, protocol_number: int, initial_count: int, final_count: int, deactivated_count: int):
        elapsed = self._stop_timer(f"pruning_{protocol_number}")
        retention = (final_count / initial_count * 100) if initial_count > 0 else 0
        
        self.info(f"\n{self.TICK} Pruning completed in {elapsed:.4f}s")
        self.info(f"Initial: {initial_count} {self.ARROW} Final: {final_count} ({retention:.1f}% retained)"
        )
        self.info(f"Deactivated: {deactivated_count}")
        self.info("")

    def skip_pruning(self, protocol_number: int): 
        self.warning(f'{self.WARNING} Pruning skipped for Protocol Step {protocol_number}.')

    # ===
    # Analysis Events
    # ===

    def pca_analysis(self, conformer_count: int, n_clusters: Optional[int], include_hydrogen: bool, output_file: str):
        self.info("")
        self.info(f"PCA Analysis:")
        self.info(f"  Conformers: {conformer_count}")
        if n_clusters:
            self.info(f"  Clusters: {n_clusters}")
        self.info(f"  Include H: {include_hydrogen}")
        self.info(f"  Output: {output_file}")
    
    def spectra_generation(self, graph_type: str, output_file: str):
        self.debug(f"  Generated {graph_type} spectrum {self.ARROW} {output_file}")

    # ===
    # Checkpoint Events
    # ===
    
    def checkpoint_saved(self, conformer_count: int):
        self.debug(f"Checkpoint saved: {conformer_count} conformers ")
    
    def checkpoint_loaded(self, conformer_count: int, protocol_number: int):
        self.info(f"Checkpoint loaded: {conformer_count} conformers")
        self.info(f"Resuming from protocol {protocol_number}")

    # ===
    # Spectra Events
    # ===

    def spectra_start(self, protocol_number:int):
        self._separator("Spectra convolution", width=40, char="-")
        self.debug(f"Starting spectra convolution for protocol {protocol_number}")
        self._start_timer(f"spectra_{protocol_number}")

    def spectra_end(self, protocol_number:int):
        elapsed = self._stop_timer(f"spectra_{protocol_number}")
        self.info(f"\n{self.TICK} Sprectra convolution completed in {elapsed:.4f}s")
        self.info("")

    def spectra_skip(self, graph_type: str): 
        self.warning(f'{self.WARNING} No calculation of {graph_type} graphs. Skipping')

    def converter_str(self, i): 
        if isinstance(i, float): 
            return f"{i:.2f}"
        return f"{i}"

    def spectra_result(self, graph_type: str, parameters: Dict, msg:str): 
        res = [f"{k}: {self.converter_str(v)}" for k, v in parameters.items()]
        self._separator(f"{graph_type} Spectra convolution", char="-", width=35)
        self.info(msg)
        self.info("\t".join(res))


    # ===
    # Error Handling
    # ===

    def critical_error(self, error_type: str, message: str, **context):
        self._separator("CRITICAL ERROR", char="!")
        self.critical(f"Error Type: {error_type}")
        self.critical(f"Message: {message}")

    # ===
    # Performance Tracking
    # ===

    def _start_timer(self, key: str):
        self._timers[key] = time.perf_counter()
    
    def _stop_timer(self, key: str) -> float:
        if key not in self._timers:
            return 0.0
        elapsed = time.perf_counter() - self._timers[key]
        del self._timers[key]
        return elapsed

    # ===
    # Formatting
    # ===

    def _separator(self, title: str = "", char: str = "=", width: int = 70, **kwargs):
        """Print a separator line."""
        if title:
            self.info(f"\n{char * width}")
            self.info(f"{title:^{width}}")
            self.info(f"{char * width}")
        else:
            self.info(f"{char * width}")
    
    def table(self, title: str, data: List[Any], headers: Union[List[str], str], **kwargs):
        """
        Log tabulated data (backward compatible with tabulate).
        
        Args:
            data: List of rows
            headers: List of column headers
        """
        self._separator(title, **kwargs)
        self.info("\n" + tabulate(data, headers=headers, floatfmt=".10f", disable_numparse=True))


    @contextmanager
    def track_operation(self, operation_name: str, **context):
        """
        Context manager to track operation duration.
        """
        start = time.perf_counter()
        context_str = ", ".join(f"{k}={v}" for k, v in context.items())
        self.debug(
            f"Starting {operation_name}" + 
            (f" ({context_str})" if context_str else "")
        )
        
        try:
            yield
        except Exception as e:
            elapsed = time.perf_counter() - start
            self.error(
                f"Operation '{operation_name}' failed after {elapsed:.2f}s: {e}"
            )
            raise
        else:
            elapsed = time.perf_counter() - start
            self.debug(f"Completed {operation_name} in {elapsed:.2f}s")
    