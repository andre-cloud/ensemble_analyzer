#!/usr/bin/env python3
import shutil
import os
import sys
import importlib.util
from pathlib import Path

# Colori per l'output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
RESET = '\033[0m'

def log_pass(msg):
    print(f"[{GREEN}PASS{RESET}] {msg}")

def log_fail(msg):
    print(f"[{RED}FAIL{RESET}] {msg}")

def log_warn(msg):
    print(f"[{YELLOW}WARN{RESET}] {msg}")

def check_python_dependencies():
    print(f"\n{'-'*20} 1. Checking Python Dependencies {'-'*20}")
    required = ['numpy', 'scipy', 'matplotlib', 'ase', 'numba', 'sklearn']
    all_pass = True
    
    for lib in required:
        if importlib.util.find_spec(lib) is not None:
            log_pass(f"Library found: {lib}")
        else:
            log_fail(f"Library MISSING: {lib}")
            all_pass = False
    
    # Check EnAn package itself
    if importlib.util.find_spec("ensemble_analyzer") is not None:
        log_pass("Ensemble Analyzer package is installed and importable.")
    else:
        log_fail("Ensemble Analyzer package NOT found. Did you run 'pip install .'?")
        all_pass = False
        
    return all_pass

def check_orca():
    print(f"\n{'-'*20} 2. Checking ORCA Configuration {'-'*20}")
    
    # Check 1: Executable in PATH
    orca_path = shutil.which("orca")
    if orca_path:
        log_pass(f"ORCA executable found at: {orca_path}")
    else:
        log_fail("ORCA executable NOT found in PATH.")
        return False

    # Check 2: Environment Variable
    version_env = os.environ.get("ORCAVERSION")
    if version_env:
        log_pass(f"ORCAVERSION environment variable is set to: {version_env}")
    else:
        log_fail("ORCAVERSION environment variable is MISSING.")
        print(f"   {YELLOW}Hint: export ORCAVERSION='x.y.z' in your shell config.{RESET}")
        return False
        
    return True

def check_gaussian():
    print(f"\n{'-'*20} 3. Checking Gaussian Configuration {'-'*20}")
    
    g16_path = shutil.which("g16")
    g09_path = shutil.which("g09")
    
    if g16_path:
        log_pass(f"Gaussian 16 found at: {g16_path}")
    elif g09_path:
        log_pass(f"Gaussian 09 found at: {g09_path}")
    else:
        log_warn("Gaussian executable (g16/g09) NOT found. (Optional if using ORCA)")

def main():
    print(f"Running Ensemble Analyzer Installation Check...\n")
    
    deps_ok = check_python_dependencies()
    orca_ok = check_orca()
    check_gaussian()
    
    print(f"\n{'-'*50}")
    if deps_ok and orca_ok:
        print(f"{GREEN}SUCCESS: Installation looks correct! You are ready to run EnAn.{RESET}")
        sys.exit(0)
    else:
        print(f"{RED}FAILURE: Critical issues found. Please fix the errors above.{RESET}")
        sys.exit(1)

if __name__ == "__main__":
    main()
