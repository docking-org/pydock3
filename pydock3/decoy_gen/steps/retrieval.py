import os
import random
import shutil
from typing import Dict, List, Tuple

from rdkit import Chem

from pydock3.decoy_gen.steps.base import DecoyGenStep
from pydock3.decoy_gen.utils import get_molecular_properties, map_to_zinc_tranche, get_progressive_windows_from_config, compare_properties_with_windows

class RetrievalStep(DecoyGenStep):
    """Retrieval step: queries ZINC20 database for property-matched decoys"""
    
    def __init__(self, working_dir: str, config):
        super().__init__(working_dir, config, "retrieval")
        self.zinc20_base_path = config.param_dict['database']['zinc20_2d_path']
        
    def run(self) -> bool:
        """Execute retrieval step"""
        try:
            self.log_info("Starting retrieval step")
            
            # Verify setup step outputs exist
            setup_dir = os.path.join(os.path.dirname(self.step_dir.path), "setup")
            ligand_map_file = os.path.join(setup_dir, "LIGAND_MAP.txt")
            if not os.path.exists(ligand_map_file):
                self.log_error("Setup step must be completed before retrieval")
                return False
            
            # Read ligand map from setup step
            ligand_map = self._read_ligand_map()
            if not ligand_map:
                self.log_error("No ligands found from setup step")
                return False
                
            self.log_info(f"Found {len(ligand_map)} ligands to process")
            
            # Process each ligand
            success_count = 0
            for ligand_num, (smiles, lig_id) in ligand_map.items():
                self.log_info(f"Processing ligand {ligand_num}: {lig_id}")
                
                success = self._query_zinc_for_ligand(ligand_num, smiles, lig_id)
                if success:
                    success_count += 1
                else:
                    self.log_error(f"Failed to retrieve decoys for ligand {ligand_num}")
            
            self.log_info(f"Successfully processed {success_count}/{len(ligand_map)} ligands")
            
            if success_count == 0:
                self.log_error("No ligands were successfully processed")
                return False
            
            self.log_info("Retrieval step completed successfully")
            return True
            
        except Exception as e:
            self.log_error(f"Retrieval step failed: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if retrieval step outputs exist"""
        # Check if we have decoy files for each ligand
        ligand_map = self._read_ligand_map()
        if not ligand_map:
            return False
            
        for ligand_num in ligand_map.keys():
            decoy_file = self.get_output_file(f"{ligand_num}_decoys.smi")
            if not os.path.exists(decoy_file):
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
    
    def _query_zinc_for_ligand(self, ligand_num: str, smiles: str, lig_id: str) -> bool:
        """Query ZINC20 database for a single ligand using progressive windowing"""
        try:
            # Calculate ligand properties
            lig_props = get_molecular_properties(smiles, get_charge=True)
            mw, logp, rotb, hbd, hba, charge = lig_props
            
            self.log_debug(f"Ligand properties - MW: {mw:.1f}, LogP: {logp:.2f}, RotB: {rotb}, HBA: {hba}, HBD: {hbd}, Charge: {charge}")
            
            # Determine ZINC20 tranche based on properties
            tranche = map_to_zinc_tranche(mw, logp)
            if not tranche:
                self.log_error(f"Could not map ligand properties to ZINC20 tranche")
                return False
                
            self.log_debug(f"Mapped to ZINC20 tranche: {tranche}")
            
            # Query ZINC20 database using progressive windowing (like original zinc_subfunc)
            decoys = self._search_zinc_tranche_progressive(tranche, lig_props)
            
            if not decoys:
                self.log_error(f"No decoys found for ligand {lig_id}")
                return False
                
            self.log_info(f"Found {len(decoys)} decoys for ligand {lig_id}")
            
            # Write decoys to file
            decoy_file = self.get_output_file(f"{ligand_num}_decoys.smi")
            with open(decoy_file, 'w') as f:
                for decoy_smiles, zinc_id, window in decoys:
                    f.write(f"{decoy_smiles} {zinc_id}\n")
                    
            return True
            
        except Exception as e:
            self.log_error(f"Error querying ZINC20 for ligand {lig_id}: {e}")
            return False
    
    def _search_zinc_tranche_progressive(self, tranche: str, lig_props: Tuple) -> List[Tuple[str, str, int]]:
        """
        Search ZINC20 tranche using progressive windowing like original zinc_subfunc
        
        Returns:
            List of (smiles, zinc_id, window_num) tuples
        """
        # Check if tranche directory exists
        tranche_dir = os.path.join(self.zinc20_base_path, tranche)
        if not os.path.exists(tranche_dir):
            self.log_debug(f"ZINC20 tranche directory not found: {tranche_dir}")
            return []
        
        # Read all SMILES files in tranche
        smi_files = [f for f in os.listdir(tranche_dir) if f.endswith('.smi')]
        if not smi_files:
            self.log_debug(f"No SMILES files found in tranche {tranche}")
            return []
        
        # Get progressive windows from config
        windows = get_progressive_windows_from_config(self.config.param_dict)
        max_decoys = self.config.param_dict['generation']['total_decoys_to_generate']
        num_windows = len(windows)
        
        # Collect all SMILES lines from all files
        all_smiles_lines = []
        for smi_file in smi_files:
            file_path = os.path.join(tranche_dir, smi_file)
            
            with open(file_path, 'r') as f:
                # Skip header if present
                first_line = f.readline().strip()
                if not first_line or first_line.lower().startswith('smiles'):
                    lines = f.readlines()
                else:
                    f.seek(0)  # Reset to beginning if no header
                    lines = f.readlines()
                all_smiles_lines.extend(lines)
        
        # Shuffle all SMILES lines using configured seed
        random.seed(self.config.param_dict['generation']['random_seed'])
        random.shuffle(all_smiles_lines)
        
        selected_decoys = []
        decoy_dict = {}
        decoy_count = 0
        best_windows = []
        
        # First pass: calculate windows for all compounds and collect window 0 matches
        self.log_debug("First pass: looking for window 0 matches...")
        for line in all_smiles_lines:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
                
            decoy_smiles, zinc_id = parts[0], parts[1]
            
            try:
                # Apply protonation filtering first
                if self._matches_protonation_state(lig_props, decoy_smiles):
                    decoy_props = get_molecular_properties(decoy_smiles, get_charge=False)
                    window = compare_properties_with_windows(lig_props, decoy_props, windows)
                    
                    if window is not None:
                        best_windows.append(window)
                        
                        # Collect window 0 matches immediately like original
                        if window == 0:
                            if zinc_id not in decoy_dict:
                                decoy_count += 1
                                decoy_dict[zinc_id] = [decoy_smiles, zinc_id, window]
                                selected_decoys.append((decoy_smiles, zinc_id, window))
                                
                                if decoy_count >= max_decoys:
                                    self.log_debug(f"DECOYS FOUND BEST WINDOW: {decoy_count}")
                                    return selected_decoys
                    else:
                        best_windows.append(None)
                else:
                    best_windows.append(None)
            except Exception as e:
                self.log_debug(f"Error processing SMILES {decoy_smiles}: {e}")
                best_windows.append(None)
                continue
        
        self.log_debug(f"DECOYS FOUND BEST WINDOW: {decoy_count}")
        
        # Second pass: try progressively wider windows like original
        for step_count in range(1, num_windows):
            self.log_debug(f"Trying window {step_count}...")
            for i, window in enumerate(best_windows):
                if window == step_count:
                    line = all_smiles_lines[i]
                    parts = line.strip().split()
                    if len(parts) < 2:
                        continue
                        
                    decoy_smiles, zinc_id = parts[0], parts[1]
                    
                    if zinc_id not in decoy_dict:
                        decoy_count += 1
                        decoy_dict[zinc_id] = [decoy_smiles, zinc_id, window]
                        selected_decoys.append((decoy_smiles, zinc_id, window))
                        
                        if decoy_count >= max_decoys:
                            self.log_debug(f"DECOYS FOUND TOTAL: {decoy_count}")
                            return selected_decoys
        
        self.log_debug(f"DECOYS FOUND TOTAL: {decoy_count}")
        return selected_decoys
    
    def run_with_scheduler(self, scheduler_name: str) -> bool:
        """Execute retrieval step using job scheduler"""
        try:
            # Import here to avoid circular imports
            from pydock3.decoy_gen.job_manager import RetrievalJobManager
            
            # Read ligand map
            ligand_map = self._read_ligand_map()
            if not ligand_map:
                self.log_error("No ligands found from setup step")
                return False
            
            self.log_info(f"Submitting {len(ligand_map)} ligand retrieval jobs to {scheduler_name}")
            
            # Submit jobs for each ligand using job manager
            job_manager = RetrievalJobManager(
                working_dir=os.path.dirname(self.step_dir.path),  # parent of retrieval dir
                config=self.config,
                scheduler_name=scheduler_name
            )
            
            success = job_manager.submit_and_wait_for_ligands(ligand_map)
            
            if success:
                self.log_info("Retrieval jobs completed successfully")
            else:
                self.log_error("All retrieval jobs failed")
                
            return success
            
        except Exception as e:
            self.log_error(f"Scheduled retrieval failed: {str(e)}")
            return False
    
    def _matches_protonation_state(self, lig_props: Tuple, decoy_smiles: str) -> bool:
        """
        Check if decoy can form a protomer that matches the ligand's charge state
        
        Args:
            lig_props: (mw, logp, rotb, hbd, hba, charge) for ligand
            decoy_smiles: SMILES string for decoy
            
        Returns:
            True if decoy can match ligand's protonation state, False otherwise
        """
        try:
            from pydock3.protonation import generate_protomers, ProtonationError
            from pydock3.decoy_gen.utils import get_progressive_windows_from_config, compare_properties_with_windows
            
            # Get ligand charge
            ligand_charge = lig_props[5]
            
            # If ligand is neutral, check if decoy is also neutral (quick check)
            if ligand_charge == 0:
                decoy_props = get_molecular_properties(decoy_smiles, get_charge=True)
                if decoy_props[5] == 0:
                    return True  # Both neutral, should match
            
            # For charged ligands or when we need to check protonation, use ChemAxon
            try:
                # Generate protomers for the decoy at pH 7.4
                smiles_with_name = [f"{decoy_smiles} temp_decoy"]
                protomer_results = generate_protomers(smiles_with_name, ph=7.4, score_cutoff=10.0)
                
                # Check if any protomer has the same charge as the ligand
                for protomer_smiles, name, score in protomer_results:
                    protomer_props = get_molecular_properties(protomer_smiles, get_charge=True)
                    protomer_charge = protomer_props[5]
                    
                    if protomer_charge == ligand_charge:
                        # Check if protomer also matches other properties
                        windows = get_progressive_windows_from_config(self.config.param_dict)
                        window = compare_properties_with_windows(lig_props, protomer_props, windows)
                        if window is not None:
                            return True
                
                return False  # No matching protomer found
                
            except (ProtonationError, ImportError):
                # If protonation tools not available, fall back to basic charge check
                decoy_props = get_molecular_properties(decoy_smiles, get_charge=True)
                return decoy_props[5] == ligand_charge
                
        except Exception as e:
            self.log_debug(f"Error in protonation matching for {decoy_smiles}: {e}")
            return False  # Conservative: reject on error


class SingleLigandRetrievalStep(DecoyGenStep):
    """Specialized step for retrieving decoys for a single ligand"""
    
    def __init__(self, working_dir: str, config, ligand_num: str, ligand_smiles: str, ligand_id: str):
        # Create unique step directory for this ligand
        step_dir = os.path.join(working_dir, "retrieval", f"job_{ligand_num}")
        super().__init__(step_dir, config, f"retrieval_{ligand_num}")
        
        self.ligand_num = ligand_num
        self.ligand_smiles = ligand_smiles
        self.ligand_id = ligand_id
        
    def run(self) -> bool:
        """Execute retrieval for this specific ligand"""
        try:
            # Use existing retrieval logic for single ligand
            # Create a temporary retrieval step instance for this ligand only
            main_working_dir = os.path.dirname(os.path.dirname(self.step_dir.path))  # Go up from retrieval/job_X to main working dir
            retrieval_step = RetrievalStep(main_working_dir, self.config)
            
            # Create a temporary ligand map with just this ligand
            temp_ligand_map = {self.ligand_num: (self.ligand_smiles, self.ligand_id)}
            
            success = retrieval_step._query_zinc_for_ligand(self.ligand_num, self.ligand_smiles, self.ligand_id)
            
            if success:
                # Copy output from main retrieval directory to job directory
                main_retrieval_dir = os.path.join(main_working_dir, "retrieval")
                source_file = os.path.join(main_retrieval_dir, f"{self.ligand_num}_decoys.smi")
                dest_file = os.path.join(self.step_dir.path, f"{self.ligand_num}_decoys.smi")
                
                if os.path.exists(source_file):
                    os.makedirs(self.step_dir.path, exist_ok=True)
                    shutil.copy2(source_file, dest_file)
                    self.log_info(f"Copied decoys file for ligand {self.ligand_num}")
                else:
                    self.log_error(f"Expected decoys file not found: {source_file}")
                    return False
                    
            return success
            
        except Exception as e:
            self.log_error(f"Single ligand retrieval failed: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if this ligand's retrieval is complete"""
        output_file = os.path.join(self.step_dir.path, f"{self.ligand_num}_decoys.smi")
        return os.path.exists(output_file)
    
