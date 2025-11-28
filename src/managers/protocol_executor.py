
from typing import List, Union
import time
import datetime
import os

import numpy as np

from src.conformer.conformer import Conformer
from src.protocol import Protocol
from src.logger.logger import Logger
from src.ioFile import save_snapshot
from src.graph import main_spectra
from src.clustering import perform_PCA, get_ensemble

from src.managers.calculation_config import CalculationConfig
from src.managers.checkpoint_manager import CheckpointManager
from src.managers.calculation_executor import CalculationExecutor

# from src.pruning import calculate_rel_energies, check_ensemble, boltzmann
from src.managers.pruning_manager import PruningManager

from src.constants import DEBUG, MIN_RETENTION_RATE





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
        self.pruning_manager = PruningManager(logger, self.config.include_H)
    
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
        
        protocol_elapsed = time.perf_counter() - protocol_start_time
        
        self.logger.info(
            f"\nTotal elapsed time for protocol {protocol.number}: "
            f"{datetime.timedelta(seconds=protocol_elapsed)}"
        )
        
        
        # Pruning
        initial_active = len([c for c in conformers if c.active])

        self.generate_report("Summary Before Pruning", conformers=conformers, protocol=protocol)
        
        self.logger.pruning_start(protocol.number, initial_active)
        
        self.pruning_manager.prune_ensemble(conformers=conformers, protocol=protocol)
        self.pruning_manager.calculate_relative_energies(conformers=conformers, temperature=self.config.temperature, protocol=protocol)
        conformers = sorted(conformers)
        
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
        if isinstance(protocol.cluster, int) and protocol.cluster > 0:
            self.logger.debug("Starting PCA" + f" {protocol.cluster=}")
            perform_PCA(
                confs=[c for c in conformers if c.active],
                ncluster=int(protocol.cluster) if isinstance(protocol.cluster, (int, float)) else None,
                fname=f"PCA_after_pruning_protocol_{protocol.number}.png",
                title=f"PCA after pruning protocol {protocol.number}",
                log=self.logger,
                include_H=self.config.include_H,
                set_=True
            )
            conformers = get_ensemble(conformers)
        
        self.generate_report("Summary After Pruning", conformers=conformers, protocol=protocol)

        self.generate_energy_report(conformers=conformers, protocol_number=protocol.number, T=self.config.temperature)


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
            active_conformers=final_active,
            deactivated=initial_active - final_active
        )

        # If retention rate is lower than 30% and TODO: setting per disabilitarlo
        retention_rate = final_active / initial_active if initial_active > 0 else 1.0

        if retention_rate < MIN_RETENTION_RATE:
            self.logger.critical(
                f'✗ Ensemble reduce more than {(1-retention_rate)*100:.1f}%. Calculation will stop.\n'
                f'\t{self.logger.WARNING} threshold: {MIN_RETENTION_RATE*100:.0f}%. Stopping.'
            )
            raise RuntimeError()
    
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
            if conf.energies.__contains__(protocol_number=str(protocol.number)):
                continue
            
            self.calculator.execute(count, conf, protocol)
            
            # Save checkpoint after each calculation
            self.checkpoint_manager.save(conformers, self.logger)
            
            count += 1

        self.checkpoint_manager.save(conformers, self.logger, log=True)



    def generate_energy_report(self, conformers: List[Conformer], protocol_number: Union[str,int], T:float):

        CONFS = [i for i in conformers if i.active]

        dE = np.array([i.energies[str(protocol_number)]["E"] for i in CONFS])
        dE_ZPVE = np.array(
            [
                i.energies[str(protocol_number)]["E"] + i.energies[str(protocol_number)]["zpve"]
                for i in CONFS
            ]
        )
        dH = np.array(
            [
                i.energies[str(protocol_number)]["E"] + i.energies[str(protocol_number)]["H"]
                for i in CONFS
            ]
        )
        dG = np.array([i.energies[str(protocol_number)]["G"] for i in CONFS])

        # Boltzmann populations
        _, dE_boltz = self.pruning_manager._boltzmann_distribution(dE, T)
        _, dE_ZPVE_boltz = self.pruning_manager._boltzmann_distribution(dE_ZPVE, T)
        _, dH_boltz = self.pruning_manager._boltzmann_distribution(dH, T)
        _, dG_boltz = self.pruning_manager._boltzmann_distribution(dG, T)

        averages = [[
            f'{T:.2f}',
            float(np.sum(dE * dE_boltz)),
            float(np.sum(dE_ZPVE * dE_ZPVE_boltz)),
            float(np.sum(dH * dH_boltz)),
            float(np.sum(dG * dG_boltz)),
        ]]

        rows = [
            [
                f"Conf {i}",
                dE[idx],
                f"{float(dE_boltz[idx]*100):.2f}",
                dE_ZPVE[idx],
                f"{float(dE_ZPVE_boltz[idx]*100):.2f}",
                dH[idx],
                f"{float(dH_boltz[idx]*100):.2f}",
                dG[idx],
                f"{float(dG_boltz[idx]*100):.2f}",
            ]
            for idx, i in enumerate(CONFS)
        ]

        self.logger.debug(rows)

        headers=["Conformer", "∆E [Eh]", "Boltzamnn Pop. on ∆E", "∆(E+ZPVE) [Eh]", "Boltzamnn Pop. on ∆(E+ZPVE)", "∆H [Eh]", "Boltzamnn Pop. on ∆H", "∆G [Eh]", "Boltzamnn Pop. on ∆G"]

        self.logger.table(
            title="Energetic Summary of the active conformers", 
            data= rows, 
            headers=headers,
            width=50, 
            char = '*'
        )

        headers=["T [K]", "E_av [Eh]", "E+ZPVE_av [Eh]", "H_av [Eh]", "G_av [Eh]"]
        self.logger.table(
            title="Ensemble Average Energies", 
            data=averages,
            headers=headers, 
            width=50, 
            char = '*',
        )

        return


    def generate_report(self, title:str, conformers: List[Conformer], protocol: Protocol):
        headers = ["Conformers",
        "E [Eh]",
        "G-E [Eh]",
        "G [Eh]",
        "B [cm-1]",
        "∆G [kcal/mol]",
        "Pop [%]",
        "Elap. time [sec]",
        "# Cluster",] + [i for i in list(protocol.verbal_internals())]

        rows = [i.create_log(protocol_number=protocol.number, monitor_internals=protocol.monitor_internals) for i in conformers if i.active]
        
        self.logger.table(title=title, data=rows, headers=headers, witdh=50, char = '*')