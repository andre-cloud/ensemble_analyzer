
from typing import List
import time
import datetime
import os

from src.conformer import Conformer
from src.protocol import Protocol
from src.logger.logger import Logger
from src.ioFile import save_snapshot
from src.pruning import calculate_rel_energies, check_ensemble
from src.graph import main_spectra, plot_comparative_graphs
from src.clustering import perform_PCA, get_ensemble

from src.managers.calculation_config import CalculationConfig
from src.managers.manager_checkpoint import CheckpointManager
from src.managers.calculation_excetuer import CalculationExecutor

from src.constants import DEBUG





class ProtocolExecutor:
    """
    Executes a single protocol on the entire ensemble.
    
    Responsibilities:
    - Run calculations for all active conformers
    - Perform pruning
    - Generate graphs and PCA
    - Save snapshots
    """
    
    def __init__(
        self,
        config: CalculationConfig,
        logger: Logger,
        checkpoint_manager: CheckpointManager
    ):
        self.config = config
        self.logger = logger
        self.checkpoint_manager = checkpoint_manager
        self.calculator = CalculationExecutor(config, logger)
    
    def execute(
        self,
        conformers: List[Conformer],
        protocol: Protocol
    ) -> None:
        """
        Execute protocol on ensemble.
        
        Args:
            conformers: Ensemble to process
            protocol: Protocol to execute
        """
        active_count = len([c for c in conformers if c.active])
        
        # Start protocol
        self.logger.protocol_start(
            number=protocol.number,
            level=protocol.calculation_level,
            functional=protocol.functional,
            basis=protocol.basis,
            active_conformers=active_count
        )
        
        protocol_start_time = time.perf_counter()
        
        # Pre-pruning PCA (if DEBUG)
        if DEBUG and protocol.cluster:
            perform_PCA(
                confs=[c for c in conformers if c.active],
                ncluster=protocol.cluster if isinstance(protocol.cluster, int) else None,
                fname=f"PCA_before_pruning_protocol_{protocol.number}.png",
                title=f"PCA before pruning protocol {protocol.number}",
                log=self.logger,
                set_=False,
                include_H=self.config.include_H
            )
        
        # Run calculations
        self._run_calculations(conformers, protocol)
        
        # Sort by energy
        conformers.sort()
        
        protocol_elapsed = time.perf_counter() - protocol_start_time
        
        self.logger.info(
            f"\nTotal elapsed time for protocol {protocol.number}: "
            f"{datetime.timedelta(seconds=protocol_elapsed)}"
        )
        
        # Pruning
        initial_active = len([c for c in conformers if c.active])
        self.logger.pruning_start(protocol.number, initial_active)
        
        check_ensemble(conformers, protocol, self.logger, self.config.include_H)
        conformers = self._sort_conformers_by_energy(conformers)
        
        final_active = len([c for c in conformers if c.active])
        
        self.logger.pruning_summary(
            protocol_number=protocol.number,
            initial_count=initial_active,
            final_count=final_active,
            deactivated_count=initial_active - final_active
        )
        
        # Save snapshot
        save_snapshot(f"ensemble_after_{protocol.number}.xyz", conformers, self.logger)
        
        # Post-pruning PCA
        if protocol.cluster:
            perform_PCA(
                confs=[c for c in conformers if c.active],
                ncluster=protocol.cluster if isinstance(protocol.cluster, int) else None,
                fname=f"PCA_after_pruning_protocol_{protocol.number}.png",
                title=f"PCA after pruning protocol {protocol.number}",
                log=self.logger,
                include_H=self.config.include_H,
                set_=True
            )
        
        # Clustering
        if isinstance(protocol.cluster, int) or protocol.cluster:
            get_ensemble(conformers)
        
        # Generate spectra
        main_spectra(
            conformers,
            protocol,
            self.logger,
            invert=self.config.invert,
            read_pop=protocol.read_population,
            fwhm=self.config.fwhm,
            shift=self.config.shift,
            definition=self.config.definition,
            interested_area=self.config.interested
        )
        
        # Protocol end
        self.logger.protocol_end(
            number=protocol.number,
            elapsed_time=datetime.timedelta(seconds=protocol_elapsed),
            active_conformers=final_active,
            deactivated=initial_active - final_active
        )
    
    def _run_calculations(
        self,
        conformers: List[Conformer],
        protocol: Protocol
    ) -> None:
        """Run calculations for all active conformers."""
        count = 1
        for conf in conformers:
            if not conf.active:
                continue
            if conf.energies.get(str(protocol.number)):
                continue
            
            self.calculator.execute(count, conf, protocol, conformers)
            
            # Save checkpoint after each calculation
            self.checkpoint_manager.save(conformers, self.logger)
            
            count += 1
    
    def _sort_conformers_by_energy(
        self,
        conformers: List[Conformer]
    ) -> List[Conformer]:
        """Sort conformers by energy."""
        calculate_rel_energies(conformers, self.config.temperature)
        return sorted(conformers)
