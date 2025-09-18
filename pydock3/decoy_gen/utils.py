"""
Utility functions for decoy generation
"""

import os
import pickle
import math
from typing import Tuple, Optional, List, Dict
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs


def get_molecular_properties(smiles: str, get_charge: bool = False) -> Tuple[float, float, int, int, int, int]:
    """
    Calculate molecular properties for a SMILES string
    
    Args:
        smiles: SMILES string
        
    Returns:
        Tuple of (molecular_weight, logp, rotatable_bonds, hb_donors, hb_acceptors, formal_charge)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    logp = Chem.Crippen.MolLogP(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    if get_charge:
        charge = Chem.GetFormalCharge(mol)
        return mw, logp, rotb, hbd, hba, charge
    else:
        return mw, logp, rotb, hbd, hba, None


def map_to_zinc_tranche(mw: float, logp: float) -> Optional[str]:
    """Map molecular weight and logP to ZINC20 tranche"""
    logp_boundaries = [
        (-1, 'A'), (0, 'B'), (1, 'C'), (2, 'D'), 
        (2.5, 'E'), (3, 'F'), (3.5, 'G'), (4, 'H'), 
        (4.5, 'I'), (5, 'J'), (float('inf'), 'K')
    ]
    
    mw_boundaries = [
        (200, 'A'), (250, 'B'), (300, 'C'), (325, 'D'), 
        (350, 'E'), (375, 'F'), (400, 'G'), (425, 'H'), 
        (450, 'I'), (500, 'J'), (float('inf'), 'K')
    ]
    
    logp_code = _find_boundary_code(logp, logp_boundaries)
    mw_code = _find_boundary_code(mw, mw_boundaries)
    
    if logp_code and mw_code:
        return f"{mw_code}{logp_code}"
    return None


def _find_boundary_code(value: float, boundaries: List[Tuple[float, str]]) -> Optional[str]:
    """Find boundary code for a given value"""
    for boundary, code in boundaries:
        if value <= boundary:
            return code
    return None

def python2_round(x, ndigits=0):
    """
    Mimic Python 2’s round(), which rounds halves away from zero
    and returns a float.
    """
    factor = 10.0 ** ndigits
    if x >= 0:
        return math.floor(x * factor + 0.5) / factor
    else:
        return math.ceil(x * factor - 0.5) / factor

def get_progressive_windows_from_config(config_dict: dict) -> List[List[float]]:
    """
    Generate progressive property windows from config ranges using linear interpolation
    
    Args:
        config_dict: Configuration dictionary containing property_matching section
        
    Returns:
        List of windows as: [mw_tol, logp_tol, rotb_tol, hbd_tol, hba_tol]
        Note: Charge matching is always exact (no tolerance)
    """
    prop_config = config_dict['property_matching']
    num_windows = prop_config['num_property_windows']
    
    # Extract min/max for each property
    mw_range = prop_config['molecular_weight_range']
    logp_range = prop_config['logp_range'] 
    rotb_range = prop_config['rotatable_bonds_range']
    hba_range = prop_config['hb_acceptors_range']
    hbd_range = prop_config['hb_donors_range']
    
    windows = []
    for i in range(num_windows):
        # Linear interpolation from min to max
        factor = i / (num_windows - 1) if num_windows > 1 else 1.0
        
        mw_tol = mw_range[0] + factor * (mw_range[1] - mw_range[0])
        logp_tol = logp_range[0] + factor * (logp_range[1] - logp_range[0])
        rotb_tol = rotb_range[0] + factor * (rotb_range[1] - rotb_range[0])
        hba_tol = hba_range[0] + factor * (hba_range[1] - hba_range[0])
        hbd_tol = hbd_range[0] + factor * (hbd_range[1] - hbd_range[0])
        
        # Ensure discrete properties remain integers
        rotb_tol = int(python2_round(rotb_tol))
        hba_tol = int(python2_round(hba_tol))
        hbd_tol = int(python2_round(hbd_tol))
        # Round logP to 1 decimal and MW to whole number
        logp_tol = python2_round(logp_tol,1)
        mw_tol = python2_round(mw_tol)

        windows.append([mw_tol, logp_tol, rotb_tol, hbd_tol, hba_tol])
    
    # Match previous functionality where MWT and logP windows won't start at 0
    if windows[0][0] == 0:
        windows[0][0] = windows[1][0]
    if windows[0][1] == 0:
        windows[0][1] = windows[1][1]

    return windows


def compare_properties_with_windows(lig_props: Tuple, dec_props: Tuple, 
                                  windows: List[List[float]]) -> Optional[int]:
    """
    Compare ligand and decoy properties using progressive property windows
    
    Args:
        lig_props: (mw, logp, rotb, hbd, hba, charge) for ligand
        dec_props: (mw, logp, rotb, hbd, hba, charge) for decoy  
        windows: Progressive tolerance windows from get_progressive_windows_from_config()
        
    Returns:
        Window number (0=best match) or None if no match
        Note: Charge matching is always exact (no tolerance)
    """
    lig_mw, lig_logp, lig_rotb, lig_hbd, lig_hba, lig_charge = lig_props
    dec_mw, dec_logp, dec_rotb, dec_hbd, dec_hba, dec_charge = dec_props
    
    # Rounding to be consistent with the old scripts
    lig_mw_round = int(python2_round(lig_mw))
    lig_logp_round = python2_round(lig_logp,2)

    for window_num, (mw_tol, logp_tol, rotb_tol, hbd_tol, hba_tol) in enumerate(windows):
        if (abs(dec_mw - lig_mw_round) <= mw_tol and
            abs(dec_logp - lig_logp_round) <= logp_tol and
            abs(dec_rotb - lig_rotb) <= rotb_tol and
            abs(dec_hbd - lig_hbd) <= hbd_tol and
            abs(dec_hba - lig_hba) <= hba_tol and
            (dec_charge is None or lig_charge is None or dec_charge == lig_charge)):
            return window_num
    
    return None

# TODO: use bulk tanimoto to make this faster
def calculate_tanimoto_matrix(ligand_smiles: List[str], decoy_smiles: List[str]) -> np.ndarray:
    """
    Calculate Tanimoto similarity matrix between ligands and decoys
    
    Args:
        ligand_smiles: List of ligand SMILES
        decoy_smiles: List of decoy SMILES
        
    Returns:
        Matrix where matrix[i][j] = Tanimoto(ligand_i, decoy_j)
    """
    # Generate fingerprints for ligands
    lig_fps = []
    for smiles in ligand_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            lig_fps.append(fp)
        else:
            lig_fps.append(None)
    
    # Generate fingerprints for decoys
    dec_fps = []
    for smiles in decoy_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            dec_fps.append(fp)
        else:
            dec_fps.append(None)
    
    # Calculate similarity matrix
    matrix = np.zeros((len(ligand_smiles), len(decoy_smiles)))
    for i, lig_fp in enumerate(lig_fps):
        for j, dec_fp in enumerate(dec_fps):
            if lig_fp is not None and dec_fp is not None:
                matrix[i][j] = DataStructs.TanimotoSimilarity(lig_fp, dec_fp)
            else:
                matrix[i][j] = 0.0
    
    return matrix


def save_tanimoto_matrix(matrix: np.ndarray, ligand_ids: List[str], decoy_ids: List[str], 
                        filepath: str) -> None:
    """Save Tanimoto matrix with metadata for caching"""
    data = {
        'matrix': matrix,
        'ligand_ids': ligand_ids,
        'decoy_ids': decoy_ids,
        'timestamp': os.path.getmtime(__file__)  # Use utils file mtime as version
    }
    
    with open(filepath, 'wb') as f:
        pickle.dump(data, f)


def load_tanimoto_matrix(filepath: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """Load cached Tanimoto matrix"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Tanimoto matrix cache not found: {filepath}")
    
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    
    return data['matrix'], data['ligand_ids'], data['decoy_ids']


def is_tanimoto_cache_valid(filepath: str, current_ligands: List[str], 
                           current_decoys: List[str]) -> bool:
    """Check if cached Tanimoto matrix is still valid"""
    try:
        _, cached_ligands, cached_decoys = load_tanimoto_matrix(filepath)
        return (set(current_ligands) == set(cached_ligands) and 
                set(current_decoys) == set(cached_decoys))
    except (FileNotFoundError, pickle.UnpicklingError):
        return False