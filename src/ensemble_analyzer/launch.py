from __future__ import annotations

from pathlib import Path
from typing import Any


def main(
    ensemble: Any = None,
    protocol: Any = None,
    *,
    restart: bool = False,
    output: str = "output.out",
    disable_color: bool = False,
    **config: Any,
) -> None:
    """Run the ensemble analysis workflow.

    Two modes
    --------
    **CLI mode** — call with no arguments; parses ``sys.argv`` via argparse
    (backward-compatible, entry point unchanged).

        >>> main()

    **Programmatic mode** — pass ``ensemble`` and / or ``protocol``.

        >>> main("ensemble.xyz", "protocol.json", temperature=350.0)

    Parameters
    ----------
    ensemble
        Path to an XYZ file, a |Conformer|, or a list of them.
    protocol
        Path to a JSON protocol file, a dict (step-number keys, same format
        as the JSON file), a |Protocol|, or a list of them.
    restart
        Resume from the last checkpoint.
    output
        Log file name.
    disable_color
        Disable ANSI colour in the log.
    **config
        All remaining keyword arguments are forwarded directly to
        |CalculationConfig|_.  See its documentation for the full list of
        accepted fields.

    .. |Conformer| replace:: ``Conformer``
    .. |Protocol| replace:: ``Protocol``
    .. |CalculationConfig| replace:: ``CalculationConfig``
    .. _CalculationConfig: ensemble_analyzer._managers.calculation_config
    """
    # ── Detect mode ────────────────────────────────────────────────────
    from_cli = ensemble is None and protocol is None

    # ── CLI mode: parse sys.argv ──────────────────────────────────────
    if from_cli:
        from ensemble_analyzer._parser_arguments import parser_arguments

        args = parser_arguments()
        ensemble = args.ensemble
        protocol = args.protocol
        restart = args.restart
        base = args.output.rsplit(".", 1)[0]
        output = f"{base}_restart.out" if restart else args.output
        disable_color = not args.disable_color

    elif restart:
        base = output.rsplit(".", 1)[0]
        output = f"{base}_restart.out"

    # ── Heavy imports (deferred) ──────────────────────────────────────
    from ensemble_analyzer.protocol.protocol import Protocol, load_protocol, sort_protocols
    from ensemble_analyzer._logger.create_log import create_logger
    from ensemble_analyzer._title import title
    from ensemble_analyzer.constants import DEBUG
    from ensemble_analyzer.ensemble_io import read_ensemble, load_workflow_data
    from ensemble_analyzer._managers.protocol_manager import ProtocolManager
    from ensemble_analyzer._managers.calculation_config import CalculationConfig
    from ensemble_analyzer._managers.calculator_orchestration import CalculationOrchestrator

    # ── Logging ────────────────────────────────────────────────────────
    log = create_logger(output_file=Path(output), debug=DEBUG, disable_color=disable_color)
    log.info(title)

    # ── Resolve inputs ─────────────────────────────────────────────────
    if restart:
        conformers, protocols_out = load_workflow_data()
        start_from = ProtocolManager().load_last_completed()
    else:
        protocols_out = _normalize_protocols(protocol)
        ProtocolManager().save(protocols_out)
        conformers = _normalize_ensemble(ensemble, log)
        start_from = 0

    protocols_out = sort_protocols(protocols_out)

    # ── Build config ───────────────────────────────────────────────────
    if from_cli:
        cfg = CalculationConfig.from_args(args, start_from)
    else:
        cfg = CalculationConfig(start_from_protocol=start_from, **config)
    cfg.restart = restart

    log.application_input_received(
        config=cfg.create_log(protocols=protocols_out, conformers=len(conformers)),
    )

    CalculationOrchestrator(
        conformers=conformers, protocols=protocols_out, config=cfg, logger=log,
    ).run()


# ---------------------------------------------------------------------------
# Normalisation helpers
# ---------------------------------------------------------------------------


def _normalize_protocols(protocol: Any) -> list:
    """Accept str | dict | Protocol | list[Protocol] → list[Protocol]."""
    from ensemble_analyzer.protocol.protocol import Protocol, load_protocol

    if isinstance(protocol, Protocol):
        return [protocol]
    if isinstance(protocol, str):
        data = load_protocol(protocol)
        return [Protocol(number=k, **data[k]) for k in data]
    if isinstance(protocol, dict):
        return [Protocol(number=k, **protocol[k]) for k in protocol]
    if isinstance(protocol, list):
        return protocol
    raise TypeError(f"protocol must be str, dict, Protocol, or list[Protocol], got {type(protocol)}")


def _normalize_ensemble(ensemble: Any, log) -> list:
    """Accept str | Conformer | list[Conformer] → list[Conformer]."""
    from ensemble_analyzer.conformer.conformer import Conformer

    if isinstance(ensemble, Conformer):
        return [ensemble]
    if isinstance(ensemble, str):
        return read_ensemble(ensemble, log)
    if isinstance(ensemble, list):
        return ensemble
    raise TypeError(f"ensemble must be str, Conformer, or list[Conformer], got {type(ensemble)}")
