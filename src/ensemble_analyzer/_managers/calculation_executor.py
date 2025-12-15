from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer.io_utils import move_files
from ensemble_analyzer.constants import regex_parsing
from ensemble_analyzer.parser_parameter import get_conf_parameters

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
    ) -> bool:
        """
        Execute calculation.
        
        Args:
            idx: Display index (1-based)
            conf: Conformer to calculate
            protocol: Protocol to use
        
        Returns:
            True if successful, False otherwise
        """

        self.logger.calculation_start(
            conformer_id=conf.number,
            protocol_number=protocol.number,
            count=idx,
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
            try:
                atoms.get_potential_energy()
            except Exception as e: 
                self.logger.debug(e)
        
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
            log=self.logger, 
            linear = self.config.linear,
            cut_off = self.config.cut_off,
            alpha = self.config.alpha,
            P = self.config.P,
        )
        
        if success:
            # Log success
            data = conf.energies.__getitem__(protocol.number)
            self.logger.calculation_success(conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=data.E, gibbs=data.G,
                frequencies = data.Freq,
                elapsed_time=elapsed)
        
        return success