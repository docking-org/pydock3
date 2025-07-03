import os
import shutil
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from pydock3.decoy_gen.steps.base import DecoyGenStep
from pydock3.decoy_gen.utils import get_molecular_properties
from pydock3.files import File

class SetupStep(DecoyGenStep):
    """Setup step: reads SMILES, optionally protonates, generates fingerprints, creates directories"""
    
    def __init__(self, working_dir: str, config, smiles_file_path: str):
        super().__init__(working_dir, config, "setup")
        self.smiles_file_path = smiles_file_path
        
    def run(self) -> bool:
        """Execute setup step"""
        try:
            self.log_info("Starting setup step")
            
            # Validate input SMILES file
            if not os.path.exists(self.smiles_file_path):
                self.log_error(f"SMILES file not found: {self.smiles_file_path}")
                return False
                
            # Read ligands from SMILES file
            ligands = self._read_ligands()
            if not ligands:
                self.log_error("No valid ligands found in SMILES file")
                return False
                
            self.log_info(f"Found {len(ligands)} ligands in input file")
            
            # Copy SMILES file to working directory
            all_ligands_path = self.get_output_file("all_ligands.smi")
            shutil.copy2(self.smiles_file_path, all_ligands_path)
            
            # Process ligands (protonate if requested)
            if self.config.param_dict['input']['protonate_ligands']:
                self.log_info("Protonation enabled - processing ligands")
                processed_ligands = self._protonate_ligands(ligands)
            else:
                self.log_info("Protonation disabled - using original ligands")
                processed_ligands = {f"{lig_id}_0": [smiles] for lig_id, smiles in ligands.items()}
            
            if not processed_ligands:
                self.log_error("No ligands after processing")
                return False
                
            self.log_info(f"Processed to {len(processed_ligands)} ligand variants")
            
            # Generate RDKit fingerprints for all ligands
            fingerprints = self._generate_fingerprints(processed_ligands)
            
            # Create individual ligand directories
            self._create_ligand_directories(processed_ligands)
            
            # Save ligand map
            self._save_ligand_map(processed_ligands)
            
            # Save fingerprints
            self._save_fingerprints(fingerprints)
            
            self.log_info("Setup step completed successfully")
            return True
            
        except Exception as e:
            self.log_error(f"Setup step failed: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if setup step outputs exist"""
        required_files = [
            "all_ligands.smi",
            "LIGAND_MAP.txt",
            "ligand_fingerprints.dat"
        ]
        
        for filename in required_files:
            if not os.path.exists(self.get_output_file(filename)):
                return False
                
        # Check if at least one ligand directory exists
        ligand_dirs = [d for d in os.listdir(self.step_dir.path) 
                      if d.startswith("ligand_") and os.path.isdir(os.path.join(self.step_dir.path, d))]
        
        return len(ligand_dirs) > 0
    
    def _read_ligands(self) -> Dict[str, str]:
        """Read ligands from SMILES file"""
        ligands = {}
        
        try:
            with open(self.smiles_file_path, 'r') as f:
                for line_num, line in enumerate(f):
                    line = line.strip()
                    if not line:
                        continue
                        
                    parts = line.split()
                    if len(parts) < 2:
                        self.log_debug(f"Skipping line {line_num+1}: insufficient columns")
                        continue
                        
                    smiles, lig_id = parts[0], parts[1]
                    
                    # Validate SMILES
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        self.log_debug(f"Skipping invalid SMILES on line {line_num+1}: {smiles}")
                        continue
                        
                    ligands[lig_id] = smiles
                    
        except Exception as e:
            self.log_error(f"Error reading SMILES file: {e}")
            return {}
            
        return ligands
    
    def _protonate_ligands(self, ligands: Dict[str, str]) -> Dict[str, List[str]]:
        """Generate protomers for ligands using ChemAxon integration"""
        from pydock3.protonation.core import generate_protomers
        
        processed = {}
        
        for lig_id, smiles in ligands.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.log_error(f"Invalid SMILES for ligand {lig_id}: {smiles}")
                continue
                
            # Generate protomers using ChemAxon
            protomers = generate_protomers(smiles, ph=7.4)
            
            if not protomers:
                self.log_error(f"Protonation failed for ligand {lig_id}: {smiles}")
                continue
                
            # Use all protomers, indexed by variant number
            for i, protomer_smiles in enumerate(protomers):
                processed[f"{lig_id}_{i}"] = [protomer_smiles]
            
            self.log_debug(f"Generated {len(protomers)} protomers for ligand {lig_id}")
            
        return processed
    
    def _generate_fingerprints(self, ligands: Dict[str, List[str]]) -> Dict[str, str]:
        """Generate RDKit fingerprints for ligands"""
        from rdkit.Chem import rdMolDescriptors
        from rdkit import DataStructs
        
        fingerprints = {}
        
        for lig_id, smiles_list in ligands.items():
            smiles = smiles_list[0]  # Use first (canonical) SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
                
            # Generate Morgan fingerprint (ECFP-like)
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp_string = DataStructs.BitVectToText(fp)
            fingerprints[lig_id] = fp_string
            
        return fingerprints
    
    def _create_ligand_directories(self, ligands: Dict[str, List[str]]) -> None:
        """Create individual directories for each ligand"""
        for i, (lig_id, smiles_list) in enumerate(ligands.items(), 1):
            ligand_dir = os.path.join(self.step_dir.path, f"ligand_{i}")
            os.makedirs(ligand_dir, exist_ok=True)
            
            # Create ligand SMILES file
            smiles_file = os.path.join(ligand_dir, f"ligand_{i}.smi")
            with open(smiles_file, 'w') as f:
                f.write(f"{smiles_list[0]} {lig_id}\n")
    
    def _save_ligand_map(self, ligands: Dict[str, List[str]]) -> None:
        """Save ligand mapping file"""
        map_file = self.get_output_file("LIGAND_MAP.txt")
        
        with open(map_file, 'w') as f:
            for i, (lig_id, smiles_list) in enumerate(ligands.items(), 1):
                f.write(f"ligand_{i} {smiles_list[0]} {lig_id}\n")
    
    def _save_fingerprints(self, fingerprints: Dict[str, str]) -> None:
        """Save fingerprints to file"""
        fp_file = self.get_output_file("ligand_fingerprints.dat")
        
        with open(fp_file, 'w') as f:
            for lig_id, fp_string in fingerprints.items():
                f.write(f"{lig_id}\t{fp_string}\n")
    
