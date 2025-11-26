from src.managers.calculation_config import CalculationConfig
from src.logger.logger import Logger

from src.conformer import Conformer
from src.protocol import Protocol
from src.IOsystem import move_files
from src.regex_parsing import regex_parsing
from src.parser_parameter import get_conf_parameters

import os
from typing import List

import time


class CalculationExecutor:
    """
    Executes single conformer calculations with retry logic.
    
    Responsibilities:
    - Run quantum calculation
    - Retry on failure with exponential backoff
    - Parse results
    - Log metrics
    """
    
    def __init__(self, config: CalculationConfig, logger: Logger):
        self.config = config
        self.logger = logger
    
    def execute(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        ensemble: List[Conformer]
    ) -> bool:
        """
        Execute calculation with retry logic.
        
        Args:
            idx: Display index (1-based)
            conf: Conformer to calculate
            protocol: Protocol to use
            ensemble: Full ensemble (for checkpoint)
        
        Returns:
            True if successful, False otherwise
        """
        success = self._single_attempt(
            idx, conf, protocol, ensemble
        )
        return success
        
                
    def _single_attempt(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        ensemble: List[Conformer],
    ) -> bool:
        """Single calculation attempt."""
        self.logger.calculation_start(
            conformer_id=conf.number,
            protocol_number=protocol.number,
            cpu=self.config.cpu,
        )
        
        # Setup calculator
        calc, label = protocol.get_calculator(cpu=self.config.cpu, conf=conf)
        atoms = conf.get_ase_atoms(calc)
        
        # Run calculation
        start_time = time.perf_counter()
        
        with self.logger.track_operation(
            "Single calculation",
            conformer_id=conf.number,
            protocol_number=protocol.number
        ):
            atoms.get_potential_energy()
        
        elapsed = time.perf_counter() - start_time
        
        # Move files
        move_files(conf, protocol, label)
        
        # Parse output
        output_file = os.path.join(
            os.getcwd(),
            conf.folder,
            f"protocol_{protocol.number}",
            f'{conf.number}_p{protocol.number}_{label}.{regex_parsing[protocol.calculator]["ext"]}'
        )
        
        # Get parameters
        success = get_conf_parameters(
            conf=conf,
            number=protocol.number,
            output=output_file,
            p=protocol,
            time=elapsed,
            temp=self.config.temperature,
            log=self.logger
        )
        
        if success:
            # Log success
            self.logger.calculation_success(conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=conf.energies[str(protocol.number)]["E"],
                elapsed_time=elapsed)
        
        return success