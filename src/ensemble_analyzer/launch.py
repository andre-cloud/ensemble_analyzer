from pathlib import Path


def main() -> None:
    from ensemble_analyzer._parser_arguments import parser_arguments

    # Parse args first — --help/--version exit here without heavy imports
    args = parser_arguments()

    # Output filename setup (light, no deps needed)
    output = args.output
    if args.restart:
        base_name = ".".join(output.split(".")[:-1])
        output = f"{base_name}_restart.out"

    # Heavy imports — only reached when a real run is requested
    from ensemble_analyzer.protocol.protocol import Protocol, load_protocol, sort_protocols
    from ensemble_analyzer._logger.create_log import create_logger
    from ensemble_analyzer._title import title
    from ensemble_analyzer.constants import DEBUG
    from ensemble_analyzer.ensemble_io import read_ensemble, load_workflow_data
    from ensemble_analyzer._managers.protocol_manager import ProtocolManager
    from ensemble_analyzer._managers.calculation_config import CalculationConfig
    from ensemble_analyzer._managers.calculator_orchestration import CalculationOrchestrator

    # Initialize logging
    log = create_logger(
        output_file=Path(output), debug=DEBUG,
        disable_color=False if not args.disable_color else True,
    )
    log.info(title)

    # Load or initialize data
    if args.restart:
        conformers, protocols = load_workflow_data()
        start_from = ProtocolManager().load_last_completed()
    else:
        protocol_data = load_protocol(args.protocol)
        protocols = [Protocol(number=idx, **protocol_data[idx]) for idx in protocol_data]
        ProtocolManager().save(protocols)
        conformers = read_ensemble(args.ensemble, log)
        start_from = 0

    # Create configuration and sort protocols
    protocols = sort_protocols(protocols)
    config = CalculationConfig.from_args(args, start_from)

    log.application_input_received(
        config=config.create_log(protocols=protocols, conformers=len(conformers)),
    )

    orchestrator = CalculationOrchestrator(
        conformers=conformers, protocols=protocols, config=config, logger=log,
    )
    orchestrator.run()
