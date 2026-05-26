from pathlib import Path
from typing import Tuple, List

from ensemble_analyzer.protocol.protocol import Protocol, load_protocol, sort_protocols
from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer._parser_arguments import parser_arguments
from ensemble_analyzer.ensemble_io import read_ensemble
from ensemble_analyzer._title import title

from ensemble_analyzer.constants import DEBUG
from ensemble_analyzer.conformer.conformer import Conformer

from ensemble_analyzer.ensemble_io import load_workflow_data
from ensemble_analyzer._managers.protocol_manager import ProtocolManager
from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._managers.calculator_orchestration import CalculationOrchestrator



def main() -> None:
    """
    Main entry point for the Ensemble Analyzer application.

    Orchestrates the full workflow:
    1. Parses command line arguments.
    2. Initializes logging.
    3. Loads or creates the Protocol and Ensemble.
    4. Configures the calculation environment.
    5. Launches the CalculationOrchestrator.

    Returns:
        None
    """
    
    # 1. Parse Arguments
    args = parser_arguments()
    
    # 2. Setup output filename
    output = args.output
    if args.restart:
        base_name = ".".join(output.split(".")[:-1])
        output = f"{base_name}_restart.out"
    
    # 3. Initialize logging
    log = create_logger(output_file=Path(output),debug=DEBUG, disable_color=False if not args.disable_color else True)
    log.info(title)
    
    # 4. Load or initialize data
    if args.restart:
        conformers, protocols = load_workflow_data()
        start_from = ProtocolManager().load_last_completed()
    else:
        # 4.1 Load protocols
        protocol_data = load_protocol(args.protocol)
        protocols = [Protocol(number=idx, **protocol_data[idx]) for idx in protocol_data]
        ProtocolManager().save(protocols)
        
        # 4.2 Load ensemble
        conformers = read_ensemble(args.ensemble, log)
        start_from = 0
    
    # 5. Create configuration and sort the protocols loaded by the number
    protocols = sort_protocols(protocols)
    config = CalculationConfig.from_args(args, start_from)
    
    # 6. Log application start
    log.application_input_received(config=config.create_log(protocols=protocols, conformers=len(conformers)))
    
    # 7. Create and run orchestrator
    orchestrator = CalculationOrchestrator(
        conformers=conformers,
        protocols=protocols,
        config=config,
        logger=log
    )
    
    # 8. Run EnAn
    orchestrator.run()
