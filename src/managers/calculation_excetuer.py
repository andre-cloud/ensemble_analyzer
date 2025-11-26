class CalculationExecutor:
    """
    Executes single conformer calculations with retry logic.
    
    Responsibilities:
    - Run quantum calculation
    - Retry on failure with exponential backoff
    - Parse results
    - Log metrics
    """
    
    def __init__(self, config: CalculationConfig, logger: EnsembleLogger):
        self.config = config
        self.logger = logger
    
    def execute(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        ensemble: List[Conformer]
    ) -> bool:
        """
        Execute calculation with retry logic.
        
        Args:
            idx: Display index (1-based)
            conf: Conformer to calculate
            protocol: Protocol to use
            ensemble: Full ensemble (for checkpoint)
        
        Returns:
            True if successful, False otherwise
        """
        for attempt in range(1, self.config.max_retries + 1):
            try:
                success = self._single_attempt(
                    idx, conf, protocol, ensemble, attempt
                )
                if success:
                    return True
                
            except Exception as e:
                self.logger.calculation_failure(
                    conformer_id=conf.number,
                    protocol_number=protocol.number,
                    error=str(e),
                    attempt=attempt,
                    will_retry=(attempt < self.config.max_retries)
                )
                
                if attempt < self.config.max_retries:
                    wait_time = min(10 * (2 ** (attempt - 1)), 300)
                    self.logger.calculation_retry(
                        conformer_id=conf.number,
                        protocol_number=protocol.number,
                        attempt=attempt + 1,
                        max_attempts=self.config.max_retries,
                        wait_seconds=wait_time
                    )
                    time.sleep(wait_time)
                else:
                    self.logger.critical_error(
                        error_type="max_retries_exceeded",
                        message=f"Max retries exceeded for CONF_{conf.number}",
                        conformer_id=conf.number,
                        protocol_number=protocol.number
                    )
                    raise RuntimeError(
                        f"Max retries ({self.config.max_retries}) exceeded"
                    ) from e
        
        return False
    
    def _single_attempt(
        self,
        idx: int,
        conf: Conformer,
        protocol: Protocol,
        ensemble: List[Conformer],
        attempt: int
    ) -> bool:
        """Single calculation attempt."""
        self.logger.calculation_start(
            conformer_id=conf.number,
            protocol_number=protocol.number,
            cpu=self.config.cpu,
            attempt=attempt
        )
        
        # Setup calculator
        calc, label = protocol.get_calculator(cpu=self.config.cpu, conf=conf)
        atoms = conf.get_ase_atoms(calc)
        
        # Run calculation
        start_time = time.perf_counter()
        
        with self.logger.track_operation(
            "quantum_calculation",
            conformer_id=conf.number,
            protocol_number=protocol.number
        ):
            atoms.get_potential_energy()
        
        elapsed = time.perf_counter() - start_time
        
        # Move files
        move_files(conf, protocol, label)
        
        # Parse output
        output_file = os.path.join(
            os.getcwd(),
            conf.folder,
            f"protocol_{protocol.number}",
            f'{conf.number}_p{protocol.number}_{label}.{regex_parsing[protocol.calculator]["ext"]}'
        )
        
        # Get parameters
        success = get_conf_parameters(
            conf=conf,
            number=protocol.number,
            output=output_file,
            p=protocol,
            time=elapsed,
            temp=self.config.temperature,
            log=self.logger
        )
        
        if success:
            # Log success
            metrics = CalculationMetrics(
                conformer_id=conf.number,
                protocol_number=protocol.number,
                energy=conf.energies[str(protocol.number)]["E"],
                elapsed_time=elapsed,
                success=True,
                attempt=attempt
            )
            self.logger.calculation_success(metrics)
        
        return success