from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import numpy as np

from collections import defaultdict

from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.protocol import Protocol
from ensemble_analyzer._logger.logger import Logger
from ensemble_analyzer.constants import R, EH_TO_KCAL, CAL_TO_J, boltzmann_distribution
from ensemble_analyzer._managers.prism_backend import PRISM_METHODS, run_prism_pruning


# ===
# Data structure
# ===

@dataclass
class ComparisonResult:
    """Result of comparing two conformers."""
    check_id: int
    reference_id: int
    delta_energy: float  # kcal/mol
    delta_rotatory: float  # cm⁻¹
    delta_moment: float  # Debye
    should_deactivate: bool
    rmsd: Optional[float] = 0
    criterion: str = "E+B"  # 'E+B' (rotconst) or 'prism:<stage> (<metric>)'
    
    def to_dict(self) -> dict:
        """
        Convert the comparison result to a dictionary for logging.

        Returns:
            dict: Dictionary representation of the comparison.
        """
        return {
            "Check": self.check_id,
            "Ref": self.reference_id,
            "∆E [kcal/mol]": f'{self.delta_energy:.2f}',
            "∆B [cm⁻¹]": f'{self.delta_rotatory:.2f}',
            "∆m [Debye]": f'{self.delta_moment:.2f}',
            "λ(D) RMSD [Å]": f'{self.rmsd:.2f}',
            "Criterion": self.criterion,
            "Deactivate": self.should_deactivate,
        }
    

class PruningManager:
    """
    Manages the ensemble refinement process (filtering and de-duplication).
    """

    def __init__(self, logger : Logger, include_H : bool = False):
        """
        Initialize the PruningManager.

        Args:
            logger (Logger): Logger instance.
            include_H (bool, optional): Whether to include Hydrogens in RMSD checks. Defaults to False.
        """

        self.logger = logger
        self.include_H = include_H
        self._deactivation_records : List[ComparisonResult] = []

    def prune_ensemble(self, conformers: List[Conformer], protocol: Protocol) -> None: 
        """
        Execute the pruning workflow on the ensemble.

        Applies:
        1. Energy window filtering (thrGMAX).
        2. Geometric de-duplication based on Energy (thrG) and Rotational Constants (thrB).

        Args:
            conformers (List[Conformer]): The full ensemble to process (modified in-place).
            protocol (Protocol): Protocol containing threshold parameters.
        """

        if self._should_skip_pruning(protocol):
            return conformers
        
        self._deactivation_records.clear() 

        # Energy window
        self._filter_by_energy_window(conformers, protocol.number, protocol.thrGMAX)

        # Structural similarity
        self._remove_duplicates(conformers, protocol)

        # Log results
        self._log_deactivations()

    def calculate_relative_energies(self, conformers: List[Conformer], temperature: float, protocol: Protocol) -> None: 
        """
        Compute relative energies and Boltzmann populations for active conformers.

        Args:
            conformers (List[Conformer]): The ensemble.
            temperature (float): Temperature [K].
            protocol (Protocol): Protocol context for energy retrieval.

        Returns:
            None
        """
        
        active = [c for c in conformers if c.active]
        if len(active) == 0: 
            return
        
        energies = np.array([c.get_energy(protocol_number=protocol.number) for c in active])
        
        rel_energies, populations = self._boltzmann_distribution(energies, temperature)

        for idx, conf in enumerate(active):
            conf.energies.last().Erel = float(rel_energies[idx])
            conf.energies.last().Pop = float(populations[idx] * 100)
        

    # ===
    # Private Functions
    # ===

    def _should_skip_pruning(self, protocol: Protocol) -> bool:
        """
        Determine whether the pruning step should be skipped.

        Pruning is skipped when the protocol is marked as no-prune.

        Args:
            protocol (Protocol): The current protocol.

        Returns:
            bool: True if pruning should be skipped, False otherwise.
        """
        if protocol.no_prune: 
            self.logger.skip_pruning(protocol_number=protocol.number)
            return True
        
        return False
    
    def _filter_by_energy_window(self, conformers: List[Conformer], protocol_number: int, threshold: float) -> None:
        """
        Deactivate conformers exceeding the maximum energy window.

        Args:
            conformers (List[Conformer]): The ensemble.
            protocol_number (int): Protocol ID to retrieve energies from.
            threshold (float): Energy window in kcal/mol.

        Returns:
            None
        """
        
        active = [(conf, self._get_effective_energy(conf)) for conf in conformers if conf.active]
        
        if len(active)==0: 
            return
        
        energies = np.array([e for _,e in active])
        rel_energies = (energies - energies.min()) * EH_TO_KCAL

        self.logger.debug(rel_energies)
        self.logger.info(f'Filtering conformers above {threshold} kcal/mol energy window')

        header = ["", "∆E [kcal/mol]"]
        rows = []
        for (conf, _), rel_e in zip(active, rel_energies):
            if rel_e > threshold:
                conf.active = False
                rows.append((f'Conf {conf.number}', f"{rel_e:.2f}"))
        
        if len(rows) > 0:
            self.logger.table(title="Conformers over energy window", data=rows, headers=header, char='*', width=30)
            self.logger.info(f"{self.logger.TICK} Deactivated {len(rows)} conformer(s)\n")
        else:
            self.logger.info(f"{self.logger.TICK} No conformers above threshold\n")


    def _remove_duplicates(self, conformers: List[Conformer], protocol: Protocol) -> None:
        """
        Identify and deactivate duplicate conformers.

        Dispatches on ``protocol.prune_method``:
        ``'rotconst'`` (default, energy + rotational constants) or
        ``'prism'`` / ``'prism_rotcorr'`` (PRISM Pruner backend).
        """
        method = getattr(protocol, "prune_method", "rotconst")
        if not isinstance(method, str):  # e.g. mocks or legacy objects
            method = "rotconst"
        if method.lower() in PRISM_METHODS:
            return self._remove_duplicates_prism(conformers, protocol, method.lower())
        return self._remove_duplicates_rotconst(conformers, protocol)

    def _remove_duplicates_prism(self, conformers: List[Conformer], protocol: Protocol, method: str) -> None:
        """
        Identify and deactivate duplicates with PRISM Pruner.

        Active conformers are sorted by effective energy, passed to
        :func:`run_prism_pruning`, and the returned mask is applied in place.
        The lowest-energy member of each duplicate pair is always retained.

        Args:
            conformers (List[Conformer]): The ensemble (modified in place).
            protocol (Protocol): Protocol with thrG and prism_* parameters.
            method (str): 'prism' or 'prism_rotcorr'.
        """
        active = [c for c in conformers if c.active]
        if len(active) < 2:
            return

        e_kcal = np.array([self._get_effective_energy(c) for c in active]) * EH_TO_KCAL
        order = np.argsort(e_kcal, kind="stable")
        active = [active[i] for i in order]
        e_kcal = e_kcal[order] - e_kcal[order].min()

        atoms = tuple(active[0].atoms)
        if any(tuple(c.atoms) != atoms for c in active):
            raise ValueError("PRISM pruning requires all conformers to share the same atom ordering")

        coords = np.array([np.asarray(c.last_geometry, dtype=float) for c in active])

        self.logger.info(f"Similarity pruning with PRISM ({method}) on {len(active)} conformers")
        keep, records = run_prism_pruning(
            coords, atoms, e_kcal,
            max_dE=float(protocol.thrG),
            method=method,
            use_moi=bool(getattr(protocol, "prism_moi", False)),
            moi_dev=float(getattr(protocol, "prism_moi_dev", 1e-2)),
            max_rmsd=float(getattr(protocol, "prism_max_rmsd", 0.25)),
            max_dev=getattr(protocol, "prism_max_dev", None),
            heavy_atoms_only=not self.include_H,
            timeout_s=int(getattr(protocol, "prism_timeout", 60)),
            debugfunction=self.logger.debug,
        )

        for rec in records:
            check = active[rec.removed]
            ref = active[rec.reference] if rec.reference >= 0 else None
            check.active = False
            check.diactivated_by = ref.number if ref is not None else None
            if ref is None:
                continue
            d_b = (abs(check.rotatory - ref.rotatory)
                   if check.rotatory is not None and ref.rotatory is not None else np.nan)
            d_m = abs(check.moment - ref.moment) if check.moment and ref.moment else np.nan
            self._deactivation_records.append(ComparisonResult(
                check_id=check.number,
                reference_id=ref.number,
                delta_energy=float(abs(e_kcal[rec.removed] - e_kcal[rec.reference])),
                delta_rotatory=d_b,
                delta_moment=d_m,
                should_deactivate=True,
                rmsd=self._calculate_rmsd(conf1=check, conf2=ref, include_H=self.include_H),
                criterion=f"prism:{rec.stage} ({rec.metric:.3f})",
            ))

        # The mask is authoritative (records may have been re-pointed by swaps).
        for conf, k in zip(active, keep):
            conf.active = bool(k)

    def _remove_duplicates_rotconst(self, conformers: List[Conformer], protocol: Protocol) -> None: 
        """
        Identify and deactivate duplicate conformers.

        Uses Rotational Constants and Energy as descriptors for fast comparison ($O(N^2)$).

        Args:
            conformers (List[Conformer]): The ensemble.
            protocol (Protocol): Protocol containing thrG and thrB thresholds.

        Returns:
            None
        """

        for idx, check in enumerate(conformers): 
            if not check.active: 
                continue
        
            for ref_idx in range(idx):
                ref = conformers[ref_idx]
                
                if ref.energies[str(protocol.number)].B == 1:
                    continue

                if not ref.active: 
                    continue

                result = self._compare_conformers(check, ref, protocol)
                # self.logger.debug(result.to_dict())
                if result.should_deactivate:
                    check.active = False
                    check.diactivated_by = ref.number
                    self._deactivation_records.append(result)
                    break

    def _compare_conformers(self, check: Conformer, ref: Conformer, protocol: Protocol) -> ComparisonResult:
        """Compare two conformers using energy and rotational constants.

        Args:
            check (Conformer): Conformer to be checked.
            ref (Conformer): Reference conformer.
            protocol (Protocol): Protocol with thresholds.

        Returns:
            ComparisonResult: Result of the comparison with deactivation decision.
        """
        delta_e = abs(self._get_effective_energy(check) - self._get_effective_energy(ref)) * EH_TO_KCAL
        delta_b = abs(check.rotatory - ref.rotatory)
        delta_m = abs(check.moment - ref.moment) if check.moment else np.nan
        
        should_deactivate = (
        delta_e < protocol.thrG and
        delta_b < protocol.thrB
        )

        comparison = ComparisonResult(
            check_id=check.number, 
            reference_id=ref.number,
            delta_energy=delta_e, 
            delta_rotatory=delta_b, 
            delta_moment=delta_m, 
            should_deactivate=should_deactivate
        )
        if should_deactivate:
            comparison.rmsd = self._calculate_rmsd(conf1=check, conf2=ref, include_H=self.include_H)
        
        return comparison
    
    # ===
    # Static Methods
    # ===

    @staticmethod
    def _get_effective_energy(conf: Conformer) -> float:
        """
        Retrieve the effective energy for a conformer.

        Args:
            conf (Conformer): The conformer to query.

        Returns:
            float: Effective energy value.
        """
        return conf.energies.get_energy()
    
    @staticmethod
    def _calculate_rmsd(conf1: Conformer, conf2: Conformer, include_H: bool) -> float:
        """Calculate RMSD based on distance matrix eigenvalues.
        It is a rotation/traslation invariant measure.

        Args:
            conf1 (Conformer): First Conformer
            conf2 (Conformer): Second Conformer
            include_H (bool): Include hydrogen atoms

        Returns:
            float: RMSD value
        """
        dm1 = conf1.distance_matrix(include_H=include_H)
        dm2 = conf2.distance_matrix(include_H=include_H)

        # Distance matrices are real-symmetric: eigvalsh returns real eigenvalues
        # sorted in ascending order, so the two spectra are compared component-wise.
        # (np.linalg.eig returns them in arbitrary order.)
        evals1 = np.linalg.eigvalsh(dm1)
        evals2 = np.linalg.eigvalsh(dm2)

        return float(np.sqrt(np.mean((evals1 - evals2) ** 2)))
    
    @staticmethod
    def _boltzmann_distribution(
        energies: np.ndarray,
        temperature: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        return boltzmann_distribution(energies, temperature)
    
    # ===
    # Logging
    # ===

    def _log_deactivations(self) -> None:
        """
        Log the deactivation records from the similarity pruning step.

        Logs a table of all conformers that were deactivated during
        duplicate removal, or a message if none were deactivated.
        """
        if not self._deactivation_records: 
            self.logger.info("No conformers deactivated by similarity check")
            return
        
        table_data = defaultdict(list)
        for record in self._deactivation_records: 
            d = record.to_dict()
            for key, value in d.items():
                table_data[key].append(value)

        self.logger.table("Conformer pruned by ∆B and ∆E", data=table_data, headers="keys", char="*", width=30)
        
        self.logger.debug(f'{self.logger.TICK} Pruning ')