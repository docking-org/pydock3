import os
from typing import Dict, List, Tuple

from pydock3.decoy_gen.steps.base import DecoyGenStep

class PreparationStep(DecoyGenStep):
    """Preparation step: generates final SMILES files for assigned decoys"""
    
    def __init__(self, working_dir: str, config):
        super().__init__(working_dir, config, "preparation")
        
    def run(self) -> bool:
        """Execute preparation step"""
        try:
            self.log_info("Starting preparation step")
            
            # Verify prerequisite step outputs exist
            filtering_dir = os.path.join(os.path.dirname(self.step_dir.path), "filtering")
            assignment_log = os.path.join(filtering_dir, "assignment_log.txt")
            
            if not os.path.exists(assignment_log):
                self.log_error("Filtering step must be completed before preparation")
                return False
            
            # Parse assignments and create SMILES files
            assignments = self._parse_assignment_log(assignment_log)
            
            if not assignments:
                self.log_error("No valid assignments found")
                return False
            
            self.log_info(f"Found assignments for {len(assignments)} ligands")
            
            # Create SMILES files for each ligand
            success_count = 0
            for lig_id, decoy_list in assignments.items():
                success = self._create_smiles_file(lig_id, decoy_list)
                if success:
                    success_count += 1
                else:
                    self.log_error(f"Failed to create SMILES file for ligand {lig_id}")
            
            # Create combined SMILES file
            combined_success = self._create_combined_smiles_file(assignments)
            
            self.log_info(f"Successfully created SMILES files for {success_count}/{len(assignments)} ligands")
            
            if success_count == 0:
                self.log_error("No SMILES files were successfully created")
                return False
            
            if not combined_success:
                self.log_error("Failed to create combined SMILES file")
                return False
            
            self.log_info("Preparation step completed successfully")
            return True
            
        except Exception as e:
            self.log_error(f"Preparation step failed: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if preparation step outputs exist"""
        # Check for combined SMILES file
        combined_file = self.get_output_file("all_decoys.smi")
        if not os.path.exists(combined_file):
            return False
        
        # Check for individual ligand SMILES files
        ligand_map = self._read_ligand_map()
        for ligand_num, (smiles, lig_id) in ligand_map.items():
            smiles_file = self.get_output_file(f"{lig_id}_decoys.smi")
            if not os.path.exists(smiles_file):
                return False
                
        return True
    
    def _read_ligand_map(self) -> Dict[str, Tuple[str, str]]:
        """Read ligand map from setup step output"""
        ligand_map = {}
        
        try:
            setup_dir = os.path.join(os.path.dirname(self.step_dir.path), "setup")
            map_file = os.path.join(setup_dir, "LIGAND_MAP.txt")
            if not os.path.exists(map_file):
                return {}
                
            with open(map_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        ligand_num = parts[0]  # e.g., "ligand_1"
                        smiles = parts[1]
                        lig_id = parts[2]
                        ligand_map[ligand_num] = (smiles, lig_id)
                        
        except Exception as e:
            self.log_error(f"Error reading ligand map: {e}")
            return {}
            
        return ligand_map
    
    def _parse_assignment_log(self, log_file: str) -> Dict[str, List[Dict]]:
        """Parse assignment log file to extract ligand-decoy assignments"""
        assignments = {}
        
        try:
            with open(log_file, 'r') as f:
                # Skip header
                header = f.readline()
                
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:  # Minimum required columns
                        lig_id = parts[0]
                        decoy_id = parts[1]
                        decoy_smiles = parts[2]
                        tc_to_lig = float(parts[3])
                        mw = float(parts[4])
                        logp = float(parts[5])
                        rotb = int(parts[6])
                        hbd = int(parts[7])
                        hba = int(parts[8]) if len(parts) > 8 else 0
                        charge = int(parts[9]) if len(parts) > 9 else 0
                        
                        if lig_id not in assignments:
                            assignments[lig_id] = []
                        
                        assignments[lig_id].append({
                            'decoy_id': decoy_id,
                            'smiles': decoy_smiles,
                            'tc_to_lig': tc_to_lig,
                            'mw': mw,
                            'logp': logp,
                            'rotb': rotb,
                            'hbd': hbd,
                            'hba': hba,
                            'charge': charge
                        })
        
        except Exception as e:
            self.log_error(f"Error parsing assignment log: {e}")
            return {}
        
        return assignments
    
    def _create_smiles_file(self, lig_id: str, decoy_list: List[Dict]) -> bool:
        """Create SMILES file for a specific ligand"""
        try:
            smiles_file = self.get_output_file(f"{lig_id}_decoys.smi")
            
            with open(smiles_file, 'w') as f:
                # Write header
                f.write("SMILES ZINC_ID TC_TO_LIG MW LogP RotB HBD HBA Charge\n")
                
                # Write decoys
                for decoy in decoy_list:
                    f.write(f"{decoy['smiles']} {decoy['decoy_id']} {decoy['tc_to_lig']:.6f} "
                           f"{decoy['mw']:.1f} {decoy['logp']:.2f} {decoy['rotb']} "
                           f"{decoy['hbd']} {decoy['hba']} {decoy['charge']}\n")
            
            self.log_info(f"Created SMILES file for {lig_id} with {len(decoy_list)} decoys")
            return True
            
        except Exception as e:
            self.log_error(f"Error creating SMILES file for {lig_id}: {e}")
            return False
    
    def _create_combined_smiles_file(self, assignments: Dict[str, List[Dict]]) -> bool:
        """Create combined SMILES file with all decoys"""
        try:
            combined_file = self.get_output_file("all_decoys.smi")
            
            with open(combined_file, 'w') as f:
                # Write header with ligand info
                f.write("SMILES ZINC_ID LIGAND_ID TC_TO_LIG MW LogP RotB HBD HBA Charge\n")
                
                total_decoys = 0
                for lig_id, decoy_list in assignments.items():
                    for decoy in decoy_list:
                        f.write(f"{decoy['smiles']} {decoy['decoy_id']} {lig_id} {decoy['tc_to_lig']:.6f} "
                               f"{decoy['mw']:.1f} {decoy['logp']:.2f} {decoy['rotb']} "
                               f"{decoy['hbd']} {decoy['hba']} {decoy['charge']}\n")
                        total_decoys += 1
            
            self.log_info(f"Created combined SMILES file with {total_decoys} total decoys from {len(assignments)} ligands")
            return True
            
        except Exception as e:
            self.log_error(f"Error creating combined SMILES file: {e}")
            return False