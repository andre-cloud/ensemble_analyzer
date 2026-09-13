import contextlib
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.protocol import Protocol
from ensemble_analyzer.ensemble_io import read_ensemble, save_snapshot
from ensemble_analyzer._managers.calculation_config import CalculationConfig
from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer.constants import DEBUG


class EnsembleResult:
    """
    Container for the results of an ensemble_analyzer workflow.
    """
    def __init__(
        self,
        conformers: List[Conformer],
        protocols: List[Protocol],
        config: CalculationConfig,
        logger: Any = None,
        orchestrator: Any = None
    ):
        self.conformers = conformers
        self.protocols = protocols
        self.config = config
        self.logger = logger
        self.orchestrator = orchestrator
        self.initial_count = len(self.conformers)

    @property
    def active(self) -> List[Conformer]:
        """List of currently active conformers."""
        return [c for c in self.conformers if c.active]

    @property
    def inactive(self) -> List[Conformer]:
        """List of inactive/pruned conformers."""
        return [c for c in self.conformers if not c.active]

    @property
    def best(self) -> Optional[Conformer]:
        """The lowest energy active conformer (assuming sorted final ensemble)."""
        active_list = self.active
        return active_list[0] if active_list else None

    @property
    def final_count(self) -> int:
        """Number of active conformers."""
        return len(self.active)

    @property
    def retention_rate(self) -> float:
        """Fraction of conformers retained."""
        return self.final_count / self.initial_count if self.initial_count > 0 else 0.0

    def __len__(self) -> int:
        return self.final_count

    def __iter__(self):
        return iter(self.active)

    def __getitem__(self, idx):
        return self.active[idx]

    def to_ase(self, active_only: bool = True) -> list:
        """Convert conformers to ASE Atoms objects."""
        confs = self.active if active_only else self.conformers
        return [c.to_ase() for c in confs]

    def save_xyz(self, path: Union[str, Path], active_only: bool = True) -> None:
        """Save the ensemble to an XYZ file."""
        confs = self.active if active_only else self.conformers
        save_snapshot(str(path), confs, log=self.logger)

    def energies(self, protocol_number: Optional[Union[int, str]] = None) -> List[Dict[str, Any]]:
        """
        Return a list of dictionaries with thermodynamic and summary data.
        
        Args:
            protocol_number: The protocol step to extract data for. If None,
                             uses the last step.
        """
        data = []
        for c in self.conformers:
            try:
                rec = c.energies.last() if protocol_number is None else c.energies[str(protocol_number)]
            except KeyError:
                rec = None
            
            data.append({
                "number": c.number,
                "active": c.active,
                "cluster": c.cluster,
                "E": getattr(rec, "E", None),
                "G": getattr(rec, "G", None),
                "H": getattr(rec, "H", None),
                "zpve": getattr(rec, "zpve", None),
                "Erel": getattr(rec, "Erel", None),
                "Pop": getattr(rec, "Pop", None),
                "B": getattr(rec, "B", None),
                "time": getattr(rec, "time", None),
            })
        return data

    def to_dataframe(self, protocol_number: Optional[Union[int, str]] = None) -> Any:
        """Convert energies data to a pandas DataFrame."""
        try:
            import pandas as pd
            return pd.DataFrame(self.energies(protocol_number=protocol_number))
        except ImportError:
            raise ImportError("pandas is required to use to_dataframe().")

    def summary(self) -> str:
        """Generate a brief textual summary of the run."""
        lines = [
            "Ensemble Analysis Summary",
            "=========================",
            f"Initial conformers: {self.initial_count}",
            f"Final active:       {self.final_count} (Retention: {self.retention_rate:.1%})",
        ]
        if self.best:
            best_e = self.best.get_energy(self.protocols[-1].number) if self.protocols else None
            lines.append(f"Lowest energy:      {best_e} (Conformer {self.best.number})")
        return "\n".join(lines)


def load_ensemble(file: Union[str, Path], raw: bool = False) -> List[Conformer]:
    """Helper to load conformers from an XYZ file without requiring a logger."""
    return read_ensemble(str(file), log=None, raw=raw)


def save_ensemble(file: Union[str, Path], conformers: List[Conformer]) -> None:
    """Helper to save conformers to an XYZ file without requiring a logger."""
    save_snapshot(str(file), conformers, log=None)


def run(
    ensemble: Any = None,
    protocol: Any = None,
    *,
    work_dir: Optional[Union[str, Path]] = None,
    restart: bool = False,
    output: Union[str, Path] = "output.out",
    quiet: bool = False,
    disable_color: bool = False,
    cpu: int = 1,
    temperature: float = 298.15,
    include_H: bool = True,
    **config: Any,
) -> EnsembleResult:
    """
    Run the complete ensemble_analyzer workflow programmatically.
    
    Args:
        ensemble: Path to XYZ file, list of Conformer, or list of ASE Atoms.
        protocol: Path to JSON protocol, dict, Protocol, list of Protocol, or None.
        work_dir: Optional directory to execute calculations within.
        restart: Resume from last checkpoint.
        output: Log file name.
        quiet: If True, suppress console logging output.
        disable_color: Disable ANSI colors in log.
        cpu: Number of CPU cores to use.
        temperature: Temperature for thermochemistry.
        include_H: Include hydrogens in analysis.
        **config: Additional CalculationConfig parameters (cut_off, alpha, linear, etc.).
        
    Returns:
        EnsembleResult: Object containing the resulting conformers and metadata.
    """
    import logging
    from ensemble_analyzer.launch import _normalize_ensemble, _normalize_protocols
    from ensemble_analyzer.protocol.protocol import sort_protocols
    from ensemble_analyzer._managers.protocol_manager import ProtocolManager
    from ensemble_analyzer._managers.calculation_config import CalculationConfig
    from ensemble_analyzer._managers.calculator_orchestration import CalculationOrchestrator
    from ensemble_analyzer.ensemble_io import load_workflow_data

    # Setup logger
    log = create_logger(
        output_file=str(output),
        debug=DEBUG,
        disable_color=disable_color
    )
    if quiet:
        # Avoid console output by raising level for handlers attached to sys.stdout
        for handler in log.handlers:
            if isinstance(handler, logging.StreamHandler):
                handler.setLevel(logging.CRITICAL)

    # Change to work_dir if requested
    if work_dir:
        Path(work_dir).mkdir(parents=True, exist_ok=True)

    cm = contextlib.chdir(work_dir) if work_dir else contextlib.nullcontext()
    with cm:
        # Load or initialize workflow
        if restart:
            conformers, protocols_out = load_workflow_data()
            start_from = ProtocolManager().load_last_completed()
        else:
            protocols_out = _normalize_protocols(protocol)
            ProtocolManager().save(protocols_out)
            conformers = _normalize_ensemble(ensemble, log)
            start_from = 0

        protocols_out = sort_protocols(protocols_out)
        
        cfg = CalculationConfig(
            cpu=cpu,
            temperature=temperature,
            include_H=include_H,
            start_from_protocol=start_from,
            **config
        )
        cfg.restart = restart

        log.application_input_received(
            config=cfg.create_log(protocols=protocols_out, conformers=len(conformers)),
        )

        orchestrator = CalculationOrchestrator(
            conformers=conformers,
            protocols=protocols_out,
            config=cfg,
            logger=log,
        )
        
        # Execute workflow (runs all protocols and finalizes)
        orchestrator.run()

        # Build and return the result
        return EnsembleResult(
            conformers=orchestrator.conformers,
            protocols=protocols_out,
            config=cfg,
            logger=log,
            orchestrator=orchestrator
        )
