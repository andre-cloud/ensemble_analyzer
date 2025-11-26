
from typing import List
import time
import datetime

from src.conformer import Conformer
from src.protocol import Protocol
from src.logger.logger import Logger
from src.ioFile import save_snapshot
from src.graph import plot_comparative_graphs
from src.clustering import perform_PCA


from src.managers.manager_checkpoint import CheckpointManager
from src.managers.protocol_manager import ProtocolManager
from src.managers.protocol_excecutor import ProtocolExecutor
from src.managers.calculation_config import CalculationConfig




class CalculationOrchestrator:
    """
    Main orchestrator for ensemble calculations.
    
    This class coordinates the entire calculation workflow:
    - Protocol execution
    - Checkpoint management
    - Final analysis and reporting
    
    This is the single entry point that replaces initiate_calculation_loop().
    """
    
    def __init__(
        self,
        conformers: List[Conformer],
        protocols: List[Protocol],
        config: CalculationConfig,
        logger: Logger
    ):
        self.conformers = conformers
        self.protocols = protocols
        self.config = config
        self.logger = logger
        
        # Managers
        self.checkpoint_manager = CheckpointManager()
        self.protocol_manager = ProtocolManager()
        
        # Executor
        self.protocol_executor = ProtocolExecutor(
            config, logger, self.checkpoint_manager
        )
    
    def run(self) -> None:
        """
        Run the complete calculation workflow.
        
        This is the main entry point that executes:
        1. Protocol loop
        2. Checkpoint management
        3. Final analysis
        4. Report generation
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
            perform_PCA(
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
        self._finalize(start_time)
    
    def _finalize(self, start_time: float) -> None:
        """Finalize calculations and generate reports."""
        # Sort final ensemble
        self.conformers = sorted(self.conformers)
        save_snapshot("final_ensemble.xyz", self.conformers, self.logger)
        
        # Final checkpoint
        self.checkpoint_manager.save(self.conformers, self.logger)
        
        # Comparative graphs
        plot_comparative_graphs(self.logger)
        
        # Calculate total time
        total_seconds = 0
        for conf in self.conformers:
            for protocol_num in conf.energies:
                total_seconds += conf.energies[protocol_num]["time"]
        
        total_time = datetime.timedelta(seconds=total_seconds)
        final_count = len([c for c in self.conformers if c.active])
        
        # Log completion
        self.logger.application_end(
            total_time=total_time,
            conformer_count=final_count
        )
        
        self.logger.info(f"\n{'='*60}")
        self.logger.info("CALCULATIONS COMPLETED SUCCESSFULLY")
        self.logger.info(f"{'='*60}\n")
