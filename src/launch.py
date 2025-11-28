from pathlib import Path
from typing import Tuple, List

from src.protocol.protocol import Protocol
from src.logger.create_log import create_logger
from src.logger.logger import Logger

from src.parser_arguments import parser_arguments
from src.protocol.protocol import load_protocol
from src.ioFile import read_ensemble
from src.title import title

from src.constants import DEBUG
from src.conformer.conformer import Conformer
from src.protocol.protocol import Protocol

from src.managers.checkpoint_manager import CheckpointManager
from src.managers.protocol_manager import ProtocolManager
from src.managers.calculation_config import CalculationConfig
from src.managers.calculator_orchestration import CalculationOrchestrator



def main():
    """
    Main entry point
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
    checkpoint_mgr = CheckpointManager()
    protocol_mgr = ProtocolManager()
    
    if args.restart:
        conformers, protocols, start_from = restart(checkpoint_mgr=checkpoint_mgr, protocol_mgr=protocol_mgr, logger=log)
    else:
        # 4.1 Load protocols
        protocol_data = load_protocol(args.protocol)
        protocols = [Protocol(number=idx, **protocol_data[idx]) for idx in protocol_data]
        protocol_mgr.save(protocols)
        
        # 4.2 Load ensemble
        conformers = read_ensemble(args.ensemble, log)
        start_from = 0
    
    # 5. Create configuration
    config = CalculationConfig.from_args(args, start_from)
    
    # 6. Log application start
    log.application_input_received({
        "temperature": config.temperature,
        "cpu": config.cpu,
        "conformers": len(conformers),
        "len_protocols": len(protocols),
        "restart": args.restart, 
        "protocols": protocols
    })
    
    # 7. Create and run orchestrator
    orchestrator = CalculationOrchestrator(
        conformers=conformers,
        protocols=protocols,
        config=config,
        logger=log
    )
    
    # 8. Run EnAn
    orchestrator.run()

def restart(
    checkpoint_mgr: CheckpointManager,
    protocol_mgr: ProtocolManager,
    logger: Logger
    ) -> Tuple[List[Conformer], List[Protocol], int]:

    conformers = checkpoint_mgr.load()
    protocols = protocol_mgr.load()
    start_from = protocol_mgr.load_last_completed()
    
    logger.checkpoint_loaded(
        conformer_count=len(conformers),
        protocol_number=start_from
    )
    
    return conformers, protocols, start_from

