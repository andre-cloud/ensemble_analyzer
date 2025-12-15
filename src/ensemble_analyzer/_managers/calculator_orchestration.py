
from typing import List
import time
import datetime

from dataclasses import dataclass

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._protocol.protocol import Protocol
from ensemble_analyzer.ensemble_io import save_snapshot
from ensemble_analyzer.graph import plot_comparative_graphs
from ensemble_analyzer.clustering import execute_PCA

from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer._managers.checkpoint_manager import CheckpointManager
from ensemble_analyzer._managers.protocol_manager import ProtocolManager
from ensemble_analyzer._managers.protocol_executor import ProtocolExecutor
from ensemble_analyzer._managers.calculation_config import CalculationConfig



@dataclass
class CalculationOrchestrator:
    """
    High-level controller for the entire analysis workflow.
    
    Manages the sequence of protocols, handles global restart logic,
    and triggers finalization steps.
    """

    conformers: List[Conformer]
    protocols: List[Protocol]
    config: CalculationConfig
    logger: Logger

    
    def __post_init__(self):        
        # Managers
        self.checkpoint_manager = CheckpointManager()
        self.protocol_manager = ProtocolManager()
        
        # Executor
        self.protocol_executor = ProtocolExecutor(self.config, self.logger, self.checkpoint_manager)
    
    def run(self) -> None:
        """
        Execute the full calculation workflow.
        
        Iterates through the protocol list starting from the configured step,
        delegating execution to ProtocolExecutor.
        """

        start_time = time.perf_counter()
        
        # Initial PCA if needed
        if len(self.conformers) > 30:
            self.logger.pca_analysis(
                conformer_count=len(self.conformers),
                n_clusters=None,
                include_hydrogen=self.config.include_H,
                output_file="initial_pca.png"
            )
            execute_PCA(
                self.conformers,
                None,
                "initial_pca.png",
                "PCA analysis of Conf Search",
                self.logger,
                set_=False,
                include_H=self.config.include_H
            )
        
        # Protocol loop
        protocols_to_run = self.protocols[self.config.start_from_protocol:]
        
        for protocol in protocols_to_run:
            # Save last protocol marker
            self.protocol_manager.save_last_completed(protocol.number)
            
            # Execute protocol
            self.protocol_executor.execute(self.conformers, protocol)
        
        # Final processing
        self._finalize(staring_time=start_time)
    
    def _finalize(self, staring_time) -> None:
        """
        Perform final wrap-up tasks after all protocols are done.

        - Saves final ensemble and checkpoint.
        - Generates comparative plots.
        - Logs total execution time and statistics.

        Args:
            staring_time (float): The workflow start time (from time.perf_counter).
        """

        # Sort final ensemble
        self.conformers = sorted(self.conformers)
        save_snapshot("final_ensemble.xyz", self.conformers, self.logger)
        
        # Final checkpoint
        self.checkpoint_manager.save(self.conformers, self.logger)
        
        # Comparative graphs
        plot_comparative_graphs(self.logger)
        
        total_time = datetime.timedelta(seconds=time.perf_counter()-staring_time)
        final_count = len([c for c in self.conformers if c.active])
        
        # Log completion
        self.logger.application_correct_end(
            total_time=total_time,
            total_conformers=final_count
        )
        
        self.logger._separator(f"CALCULATIONS COMPLETED SUCCESSFULLY", char="*")
