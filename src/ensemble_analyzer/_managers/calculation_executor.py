from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._protocol.protocol import Protocol

from ensemble_analyzer.constants import regex_parsing
from ensemble_analyzer.parser_parameter import get_conf_parameters
from ensemble_analyzer._calculators.base import ML_CALCULATORS
from ensemble_analyzer._conformer.energy_data import EnergyRecord, compute_rotational_constants

import os
import shutil
from pathlib import Path
import time
import numpy as np


class CalculationExecutor:
    """
    Executes single conformer calculations.
    
    Orchestrates the lifecycle of a single QM job: input generation,
    execution, file management, and result parsing.
    """
    
    def __init__(self, config: CalculationConfig, logger: Logger) -> None:
        """
        Initialize the executor.

        Args:
            config (CalculationConfig): Global configuration.
            logger (Logger): Application logger.
        """

        self.config = config
        self.logger = logger
    
    def execute(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        cpu: int | None = None,
    ) -> bool:
        """
        Run a calculation for a specific conformer and protocol.

        Args:
            idx (int): Display index (1-based count for logging).
            conf (Conformer): The conformer to calculate.
            protocol (Protocol): The computational protocol to apply.
            cpu (int | None): CPUs for this job. Defaults to ``self.config.cpu``.

        Returns:
            bool: True if the calculation and parsing were successful, False otherwise.
        """

        self.logger.calculation_start(
            conformer_id=conf.number,
            protocol_number=protocol.number,
            count=idx,
        )

        per_job_cpu = cpu if cpu is not None else self.config.cpu
        
        is_ml = protocol.calculator.lower() in ML_CALCULATORS
        
        # Setup calculator
        calc, label = protocol.get_calculator(cpu=per_job_cpu, conf=conf)
        atoms = conf.get_ase_atoms(calc)
        
        # Ensure output directory exists (ASE writes files via label path)
        os.makedirs(f"{conf.folder}/protocol_{protocol.number}", exist_ok=True)
        
        # Run calculation
        os.environ['OMP_NUM_THREADS'] = str(per_job_cpu)
        os.environ['MKL_NUM_THREADS'] = str(per_job_cpu) 
        os.environ['OPENBLAS_NUM_THREADS'] = str(per_job_cpu)

        start_time = time.perf_counter()

        with self.logger.track_operation(
            "Single calculation",
            conformer_id=conf.number,
            protocol_number=protocol.number
        ):
            try:
                energy = atoms.get_potential_energy()
            except Exception as e:
                self.logger.debug(e)
                return False

        elapsed = time.perf_counter() - start_time

        if is_ml:
            conf.energies.add(
                protocol.number,
                EnergyRecord(E=energy, time=elapsed),
            )
            compute_rotational_constants(conf, protocol.number)
            self.logger.calculation_success(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=energy, gibbs=np.nan,
                frequencies=np.array([]),
                elapsed_time=elapsed,
            )
            return True

        calc_name = protocol.calculator.lower()
        ext = regex_parsing[calc_name]["ext"]
        output_file = os.path.join(
            os.getcwd(),
            conf.folder,
            f"protocol_{protocol.number}",
            f'{conf.number}_p{protocol.number}_{label}.{ext}'
        )

        # GenericFileIOCalculator (ORCA, etc.) writes with template-defined
        # names inside calc.directory; rename to match parser expectation
        if hasattr(calc, 'template') and hasattr(calc.template, 'outputname'):
            src = Path(calc.directory).resolve() / calc.template.outputname
            dst = Path(output_file)
            if src.exists() and src != dst:
                shutil.move(str(src), str(dst))

        success = get_conf_parameters(
            conf=conf,
            number=protocol.number,
            output=output_file,
            p=protocol,
            time=elapsed,
            temp=self.config.temperature,
            log=self.logger,
            linear=self.config.linear,
            cut_off=self.config.cut_off,
            alpha=self.config.alpha,
            P=self.config.P,
        )

        if success:
            data = conf.energies[protocol.number]
            self.logger.calculation_success(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=data.E, gibbs=data.G,
                frequencies=data.Freq,
                elapsed_time=elapsed,
            )

        return success