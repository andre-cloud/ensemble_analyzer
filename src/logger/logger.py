import logging
import sys
import time
from pathlib import Path
from typing import Optional, Any, Dict
from contextlib import contextmanager
from datetime import timedelta
from functools import wraps

from src.constants import DEBUG, ordinal
from src.title import title



level_default = logging.DEBUG if DEBUG else logging.INFO



class Logger(logging.Logger): 

    def __init__(self, name:str, level=level_default):
        super().__init__(name, level=level)

        self._timers: Dict[str, float] = {}

    # ===
    # Application Event
    # ===
    
    def title_screen(self):
        self.info(title)

    def application_input_recieved(self, config: Dict[str,Any]): 
        self._separator("Calculation Input")
        self.info(f"Ensemble: {config.get('conformers', 'N/A')} confromer(s)")
        self.indo(f'Protocols: {config.get('protocols', 'N/A')}')
        self.info(f"Temperature: {config.get('temperature', 'N/A')} K")
        self.info(f"CPU cores: {config.get('cpu', 'N/A')}")
        if config.get('restart'):
            self.info("Mode: RESTART")
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
        self.info(f"Active: {active_conformers} | Deactivated: {deactivated}")
        self._separator(f"END PROTOCOL {number}")
        self.info("")

    # === 
    # Calculation Events
    # ===

    def calculation_start(self, conformer_id: int, protocol_number: int, cpu: int):
        self.info(f"  → CONF {conformer_id:03d} | Protocol {protocol_number} | CPU {cpu}")
        self._start_timer(f"calc_{conformer_id}_{protocol_number}")
    
    def calculation_success(self, conformer_id: int, protocol_number: int, energy: float, elapsed_time: float):
        self._stop_timer(f"calc_{conformer_id}_{protocol_number}")
        self.info(f"    ✓ CONF {conformer_id:03d} | E = {energy:.8f} Eh | Time: {elapsed_time:.1f}s")
    
    def calculation_failure(self, conformer_id: int, error: str):
        status = "FAILED"
        self.error(f"    ✗ CONF {conformer_id:03d} [{status}] | Error: \n{error[:60]}")


    # ===
    # Pruning Events
    # ===

    def pruning_start(self, protocol_number: int, conformer_count: int):
        self.info("")
        self.info(f"Starting pruning for protocol {protocol_number}")
        self.info(f"Conformers before pruning: {conformer_count}")
        self._start_timer(f"pruning_{protocol_number}")

    def conformer_deactivated(self, conformer_id: int, reason: str, reference_id: Optional[int] = None, delta_energy: Optional[float] = None, delta_rotatory: Optional[float] = None):
        ref_str = f" (ref: CONF {reference_id})" if reference_id else ""
        details = []
        if delta_energy is not None:
            details.append(f"ΔE={delta_energy:.3f} kcal/mol")
        if delta_rotatory is not None:
            details.append(f"ΔB={delta_rotatory:.3f} cm⁻¹")
        detail_str = f" [{', '.join(details)}]" if details else ""
        
        self.debug(f"  ⊗ CONF {conformer_id:03d} deactivated | {reason}{ref_str}{detail_str}")

    def pruning_summary(self, protocol_number: int, initial_count: int, final_count: int, deactivated_count: int):
        elapsed = self._stop_timer(f"pruning_{protocol_number}")
        retention = (final_count / initial_count * 100) if initial_count > 0 else 0
        
        self.info(f"Pruning completed in {elapsed:.2f}s")
        self.info(f"Initial: {initial_count} → Final: {final_count} ({retention:.1f}% retained)"
        )
        self.info(f"Deactivated: {deactivated_count}")
        self.info("")

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
        self.debug(f"  Generated {graph_type} spectrum → {output_file}")

    # ===
    # Checkpoint Events
    # ===
    
    def checkpoint_saved(self, conformer_count: int):
        self.debug(f"Checkpoint saved: {conformer_count} conformers ")
    
    def checkpoint_loaded(self, conformer_count: int, protocol_number: int):
        self.info(f"Checkpoint loaded: {conformer_count} conformers")
        self.info(f"Resuming from protocol {protocol_number}")

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

    def _separator(self, title: str = "", char: str = "=", width: int = 70):
        """Print a separator line."""
        if title:
            self.info(f"\n{char * width}")
            self.info(f"{title:^{width}}")
            self.info(f"{char * width}")
        else:
            self.info(f"{char * width}")
    
    def table(self, data, headers):
        """
        Log tabulated data (backward compatible with tabulate).
        
        Args:
            data: List of rows
            headers: List of column headers
        """
        from tabulate import tabulate
        self.info("\n" + tabulate(data, headers=headers, floatfmt=".5f"))
