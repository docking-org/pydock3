"""
Core protonation functionality for PyDock3

This module provides Python-based protomer generation using ChemAxon tools
and filtering based on protonation state matching.
"""

import subprocess
import tempfile
import os
import logging
from typing import List, Tuple, Optional, Dict
from collections import defaultdict

from pydock3.decoy_gen.utils import get_molecular_properties, compare_properties_with_windows

logger = logging.getLogger(__name__)

class ProtonationError(Exception):
    """Raised when protonation operations fail"""
    pass

def check_chemaxon_tools() -> bool:
    """
    Check if ChemAxon cxcalc and molconvert tools are available
    
    Returns:
        True if both tools are available, False otherwise
    """
    from .config import get_chemaxon_paths, setup_chemaxon_environment
    
    # Get configured paths
    cxcalc_path, molconvert_path = get_chemaxon_paths()
    
    if not (cxcalc_path and molconvert_path):
        logger.debug("ChemAxon tools not configured. Run 'pydock3 configure' to set up.")
        return False
    
    # Setup environment
    if not setup_chemaxon_environment():
        logger.debug("ChemAxon license not configured")
        return False
    
    try:
        # Test both tools
        subprocess.run([cxcalc_path, "--help"], capture_output=True, check=True, timeout=10)
        subprocess.run([molconvert_path, "-h"], capture_output=True, check=True, timeout=10)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.debug(f"ChemAxon tools check failed: {e}")
        return False

def generate_protomers(smiles_with_names: List[str], ph: float = 7.4,
                      tautomer_limit: int = 20, protomer_limit: int = 20,
                      score_cutoff: float = 10.0, timeout: int = 300) -> List[Tuple[str, str, float]]:
    """
    Generate protomers for a list of SMILES using ChemAxon tools
    
    Args:
        smiles_with_names: List of "SMILES name" strings or just "SMILES"
        ph: pH for protonation (default 7.4)
        tautomer_limit: Minimum tautomer score (default 20)
        protomer_limit: Minimum protomer score (default 20)  
        score_cutoff: Minimum combined score (default 10.0)
        timeout: Timeout in seconds for each subprocess call
        
    Returns:
        List of (protomer_smiles, original_name, combined_score)
        
    Raises:
        ProtonationError: If ChemAxon tools fail or are not available
    """
    if not check_chemaxon_tools():
        raise ProtonationError("ChemAxon cxcalc and molconvert tools are required but not available. "
                             "Run 'python -m pydock3.scripts configure' to set up ChemAxon tools.")
    
    # Get configured tool paths
    from .config import get_chemaxon_paths, setup_chemaxon_environment
    cxcalc_path, molconvert_path = get_chemaxon_paths()
    setup_chemaxon_environment()
    
    if not smiles_with_names:
        return []
    
    logger.info(f"Generating protomers for {len(smiles_with_names)} molecules at pH {ph}")
    
    # Prepare input - ensure each line has SMILES + name format
    input_lines = []
    for i, line in enumerate(smiles_with_names):
        parts = line.strip().split()
        if len(parts) >= 2:
            input_lines.append(line.strip())
        elif len(parts) == 1:
            input_lines.append(f"{parts[0]} mol_{i+1}")
        else:
            logger.warning(f"Skipping invalid SMILES line: {line}")
    
    if not input_lines:
        logger.warning("No valid SMILES found")
        return []
    
    input_text = "\n".join(input_lines)
    logger.debug(f"Input text: {input_text}")
    
    try:
        # Step 1: Generate tautomers
        logger.debug("Step 1: Generating tautomers")
        cmd1 = [cxcalc_path, "-g", "dominanttautomerdistribution", "-H", str(ph), "-C", "false", "-t", "tautomer-dist"]
        
        result1 = subprocess.run(cmd1, input=input_text, text=True, capture_output=True, timeout=timeout)
        if result1.returncode != 0:
            raise ProtonationError(f"Tautomer generation failed: {result1.stderr}")
        
        # Step 2: Filter tautomers by score
        logger.debug("Step 2: Filtering tautomers")
        cmd2 = [molconvert_path, "sdf", "-g", "-c", f"tautomer-dist>={tautomer_limit}"]
        
        result2 = subprocess.run(cmd2, input=result1.stdout, text=True, capture_output=True, timeout=timeout)
        if result2.returncode != 0:
            raise ProtonationError(f"Tautomer filtering failed: {result2.stderr}")
        
        # Step 3: Generate protomers
        logger.debug("Step 3: Generating protomers")
        cmd3 = [cxcalc_path, "-g", "microspeciesdistribution", "-H", str(ph), "-t", "protomer-dist"]
        
        result3 = subprocess.run(cmd3, input=result2.stdout, text=True, capture_output=True, timeout=timeout)
        if result3.returncode != 0:
            raise ProtonationError(f"Protomer generation failed: {result3.stderr}")
        
        # Step 4: Convert to SMILES with scores
        logger.debug("Step 4: Converting to SMILES")
        cmd4 = [molconvert_path, "smiles", "-g", "-c", f"protomer-dist>={protomer_limit}", 
                "-T", "name:tautomer-dist:protomer-dist"]
        
        result4 = subprocess.run(cmd4, input=result3.stdout, text=True, capture_output=True, timeout=timeout)
        if result4.returncode != 0:
            raise ProtonationError(f"SMILES conversion failed: {result4.stderr}")
        
        # Step 5: Parse results and filter by score
        results = []
        lines = result4.stdout.strip().split('\n')
        
        for line in lines:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    protomer_smiles = parts[0]
                    original_name = parts[1] if parts[1] else f"mol_{len(results)+1}"
                    taut_score = float(parts[2])
                    prot_score = float(parts[3])
                    
                    # Calculate combined score: (taut * prot) / 100
                    combined_score = (taut_score * prot_score) / 100.0
                    
                    if combined_score >= score_cutoff:
                        results.append((protomer_smiles, original_name, combined_score))
        
        logger.info(f"Generated {len(results)} protomers above score {score_cutoff}")
        return results
        
    except subprocess.TimeoutExpired as e:
        raise ProtonationError(f"Protomer generation timed out after {timeout}s")
    except Exception as e:
        raise ProtonationError(f"Protomer generation failed: {e}")


def protonate_smiles_file(input_file: str, output_file: str, ph: float = 7.4, 
                         score_cutoff: float = 10.0) -> None:
    """
    Protonate SMILES from a file and write results to output file
    
    Args:
        input_file: Path to input file with SMILES (format: "SMILES name" per line)
        output_file: Path to output file for results
        ph: pH for protonation
        score_cutoff: Minimum combined score
        
    Raises:
        ProtonationError: If file I/O or protomer generation fails
    """
    try:
        # Read input file
        with open(input_file, 'r') as f:
            smiles_lines = [line.strip() for line in f if line.strip()]
        
        if not smiles_lines:
            raise ProtonationError(f"No SMILES found in {input_file}")
        
        logger.info(f"Read {len(smiles_lines)} SMILES from {input_file}")
        
        # Generate protomers
        results = generate_protomers(smiles_lines, ph=ph, score_cutoff=score_cutoff)
        
        # Write output
        with open(output_file, 'w') as f:
            f.write("SMILES\tName\tScore\n")
            for protomer_smiles, name, score in results:
                f.write(f"{protomer_smiles}\t{name}\t{score:.2f}\n")
        
        logger.info(f"Wrote {len(results)} protomers to {output_file}")
        
    except IOError as e:
        raise ProtonationError(f"File I/O error: {e}")
    except Exception as e:
        raise ProtonationError(f"Protonation failed: {e}")