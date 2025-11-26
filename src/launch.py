from pathlib import Path
import json
import os

from src.protocol import Protocol
from src.logger.create_log import create_logger

from src.parser_arguments import parser_arguments
from src.protocol import load_protocol
from src.ioFile import read_ensemble
from src.title import title

from src.constants import DEBUG


from src.managers.manager_checkpoint import CheckpointManager
from src.managers.protocol_manager import ProtocolManager
from src.managers.calculation_config import CalculationConfig
from src.managers.calculator_orchestration import CalculationOrchestrator



def main():
    """
    Main entry point
    """
    
    args = parser_arguments()
    
    # Load or create settings
    if os.path.exists("settings.json"):
        settings = json.load(open("settings.json"))
    else:
        settings = {
            "output": args.output,
            "cpu": args.cpu,
            "temperature": args.temperature,
            "definition": args.definition,
            "fwhm_vibro": args.fwhm_vibro,
            "fwhm_electro": args.fwhm_electro,
            "shift_vibro": args.shift_vibro,
            "shift_electro": args.shift_electro,
            "interested_vibro": args.interest_vibro ,
            "interested_electro": args.interest_electro ,
            "invert": args.invert,
            "include_H": args.exclude_H,
        }
        json.dump(settings, open("settings.json", "w"), indent=4)
    
    # Setup output file
    output = (
        settings.get("output", args.output)
        if not args.restart
        else ".".join(settings.get("output", args.output).split(".")[:-1])
        + "_restart.out"
    )
    
    # Initialize structured logging
    logger = create_logger(output_file=Path(output),debug=DEBUG, also_console=False)
    
    logger.info(title)
    
    # Load or initialize data
    checkpoint_mgr = CheckpointManager()
    protocol_mgr = ProtocolManager()
    
    if args.restart:
        conformers = checkpoint_mgr.load()
        protocols = protocol_mgr.load()
        start_from = protocol_mgr.load_last_completed()
        
        logger.checkpoint_loaded(
            conformer_count=len(conformers),
            protocol_number=start_from
        )
    else:
        # Load protocols
        protocol_data = load_protocol(args.protocol)
        protocols = [Protocol(number=idx, **protocol_data[idx]) for idx in protocol_data]
        protocol_mgr.save(protocols)
        
        # Load ensemble
        conformers = read_ensemble(args.ensemble, logger)
        start_from = 0
    
    # Create configuration
    config = CalculationConfig(
        cpu=settings.get("cpu", args.cpu),
        temperature=settings.get("temperature", args.temperature),
        start_from_protocol=start_from,
        include_H=settings.get("include_H", True),
        definition=settings.get("definition", 4),
        fwhm={
            'vibro': settings.get("fwhm_vibro"),
            'electro': settings.get("fwhm_electro")
        },
        shift={
            'vibro': settings.get("shift_vibro"),
            'electro': settings.get("shift_electro")
        },
        interested={
            'vibro': settings.get("interested_vibro"),
            'electro': settings.get("interested_electro")
        },
        invert=settings.get("invert", False)
    )
    
    # Log application start
    logger.application_start({
        "temperature": config.temperature,
        "cpu": config.cpu,
        "conformers": len(conformers),
        "protocols": len(protocols),
        "restart": args.restart
    })
    
    # Create and run orchestrator
    orchestrator = CalculationOrchestrator(
        conformers=conformers,
        protocols=protocols,
        config=config,
        logger=logger
    )
    
    orchestrator.run()