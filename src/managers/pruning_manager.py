from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import numpy as np

from collections import defaultdict

from src.conformer import Conformer
from src.protocol import Protocol
from src.logger.logger import Logger
from src.constants import R, EH_TO_KCAL, CAL_TO_J


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
    rmsd: Optional[float] = None
    
    def to_dict(self) -> Dict:
        return {
            "Check": self.check_id,
            "Ref": self.reference_id,
            "∆E [kcal/mol]": self.delta_energy,
            "∆B [cm⁻¹]": self.delta_rotatory,
            "∆m [Debye]": self.delta_moment,
            "λi RMSD": self.rmsd,
            "Deactivate": self.should_deactivate,
        }
    

class PruningManager:

    def __init__(self, logger : Logger, include_H : bool = True): 
        self.logger = logger
        self.include_H = include_H
        self._deactivation_records : List[ComparisonResult] = []

    def prune_ensemble(self, conformers: List[Conformer], protocol: Protocol) -> List[Conformer]: 
        """Main pruning workflow

        Args:
            conformers (List[Conformer]): Ensemble to prune
            protocol (Protocol): Protocol

        Returns:
            List[Conformer]: Pruned ensemble (same list, but modified in-place)
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

    def calculate_relative_energies(self, conformers: List[Conformer], temperature: float) -> None: 
        """Calculate relative energies and Boltzmann populations.

        Args:
            conformers (List[Conformer]): Ensemble
            temperature (float): Temperature [K]
        """
        
        active = [c for c in conformers if c.active]
        if len(active) == 0: 
            return
        
        energies = np.array([c.get_energy for c in active])
        
        rel_energies, populations = self._boltzmann_distribution(energies, temperature)

        for idx, conf in enumerate(active):
            conf._last_energy["Erel"] = float(rel_energies[idx])
            conf._last_energy["Pop"] = float(populations[idx] * 100)
        

    # ===
    # Private Functions
    # ===

    def _should_skip_pruning(self, protocol: Protocol) -> bool : 
        """Check id pruning should be skipped: protocol.no_prune or protocol.graph"""
        if protocol.graph or protocol.no_prune: 
            self.logger.skip_pruning()
            return True
        
        return False
    
    def _filter_by_energy_window(self, conformers: List[Conformer], protocol_number: int, threshold: float):
        """Deactivate conformers above energy window

        Args:
            conformers (List[Conformer]): Ensemble
            protocol_number (int): Current Protocol
            threshold (float): Max energy window [kcal/mol]
        """
        
        active = [(conf, self._get_effective_energy(conf, protocol_number)) for conf in conformers if conf.active]
        
        if len(active)==0: 
            return
        
        energies = np.array([e for _,e in active])
        rel_energies = (energies - energies.min()) * EH_TO_KCAL

        self.logger.info(f'Filtering conformers above {threshold} kcal/mol energy window')

        header = ["", "∆E [kcal/mol]"]
        rows = []
        for (conf, _), rel_e in zip(active, rel_energies):
            if rel_e > threshold:
                conf.active = False
                rows.append((active, f"{rel_e:.2f}"))
        
        if len(rows) > 0:
            self.logger.table(title="Conformers over energy window", data=rows, headers=header, char='*', width=30)
            self.logger.info(f"{self.logger.TICK} Deactivated {len(rows)} conformer(s)\n")
        else:
            self.logger.info(f"{self.logger.TICK} No conformers above threshold\n")


    def _remove_duplicates(self, conformers: List[Conformer], protocol: Protocol) -> None: 
        """Remove duplicate conformers based on energy and magnitude of inertia momentum

        Args:
            conformers (List[Conformer]): Ensemble
            protocol (Protocol): Protocol with thresholds
        """

        for idx, check in enumerate(conformers): 
            if not check.active: 
                continue
        
            for ref_idx in range(idx):
                ref = conformers[ref_idx]
                
                if not ref.active: 
                    continue

                result = self._compare_conformers(check, ref, protocol)
                if result.should_deactivate:
                    check.active = False
                    check.diactivated_by = ref.number
                    self._deactivation_records.append(result)
                    break

    def _compare_conformers(self, check: Conformer, ref: Conformer, protocol: Protocol) -> ComparisonResult:
        """Compare two conformers

        Args:
            check (Conformer): Conformer to be check
            ref (Conformer): Reference conformer
            protocol (Protocol): Protocol with thresholds
        """
        delta_e = (self._get_effective_energy(check, protocol_number=protocol.number) - self._get_effective_energy(ref, protocol_number=protocol.number)) * EH_TO_KCAL
        delta_b = abs(check.rotatory - ref.rotatory)
        delta_m = abs(check.moment - ref.moment)
        
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
    def _get_effective_energy(conf: Conformer, protocol_number:str) -> float: 
        """Get G if available, otherwise E."""
        energy_data = conf.energies.get(str(protocol_number))
        if not energy_data:
            return 0.0
        
        g = energy_data.get("G")
        if g is not None and not np.isnan(g):
            return g
        return energy_data.get("E", 0.0)
    
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

        evals1, _ = np.linalg.eig(dm1)
        evals2, _ = np.linalg.eig(dm2)

        return float(np.sqrt(np.mean((evals1 - evals2) ** 2)))
    
    @staticmethod
    def _boltzmann_distribution(
        energies: np.ndarray, 
        temperature: float
    ) -> Tuple[np.ndarray, np.ndarray]: 
        """Calculate the Boltzmann distribution

        Args:
            energies (np.ndarray): Array of energies [Eh]
            temperature (float): Temperature [K]

        Returns:
            Tuple[np.ndarray, np.ndarray]: Relative_energies and population
        """

        rel_energies = energies - energies.min()
        exponent = -(rel_energies * EH_TO_KCAL * 1000 * CAL_TO_J) / (R * temperature)
        boltz_weights = np.exp(exponent)

        population = boltz_weights / boltz_weights.sum()

        return rel_energies, population
    
    # ===
    # Logging
    # ===

    def _log_deactivations(self) -> None:
        self.logger.info(f'{self.logger.TICK} Pruning ')
        if not self._deactivation_records: 
            self.logger.info("No conformers deactivated by similarity check")
            return
        
        table_data = defaultdict(list)
        for record in self._deactivation_records: 
            d = record.to_dict()
            for key, value in d.items():
                table_data[key].append(value)

        self.logger.table("Conformer pruned by ∆B and ∆E", data=table_data, headers="keys", char="*", width=30)