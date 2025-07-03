import os
import sys
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

import pulp
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

from pydock3.decoy_gen.steps.base import DecoyGenStep
from pydock3.decoy_gen.utils import (get_molecular_properties, calculate_tanimoto_matrix, 
                                    save_tanimoto_matrix, load_tanimoto_matrix, is_tanimoto_cache_valid,
                                    get_progressive_windows_from_config, compare_properties_with_windows)

class FilteringStep(DecoyGenStep):
    """Filtering step: calculates Tanimoto coefficients, performs clustering, and optimizes decoy assignment"""
    
    def __init__(self, working_dir: str, config):
        super().__init__(working_dir, config, "filtering")
        
    def run(self) -> bool:
        """Execute filtering step"""
        try:
            self.log_info("Starting filtering step")
            
            # Verify prerequisite step outputs exist
            setup_dir = os.path.join(os.path.dirname(self.step_dir.path), "setup")
            retrieval_dir = os.path.join(os.path.dirname(self.step_dir.path), "retrieval")
            
            ligand_map_file = os.path.join(setup_dir, "LIGAND_MAP.txt")
            if not os.path.exists(ligand_map_file):
                self.log_error("Setup step must be completed before filtering")
                return False
                
            # Check that retrieval step has some decoy files
            if not os.path.exists(retrieval_dir):
                self.log_error("Retrieval step must be completed before filtering")
                return False
            
            # Collect ligands and decoys with properties
            lig_property_dict, decoy_property_dict = self._collect_ligands_and_decoys()
            
            if not lig_property_dict:
                self.log_error("No ligands found")
                return False
                
            if not decoy_property_dict:
                self.log_error("No decoys found")
                return False
                
            self.log_info(f"Found {len(lig_property_dict)} ligands and {len(decoy_property_dict)} decoys")
            
            # Calculate Tanimoto coefficients (protonation filtering now done in retrieval step)
            self.log_info("Calculating Tanimoto coefficients...")
            decoy_tc_list = self._calculate_tanimoto_coefficients(lig_property_dict, decoy_property_dict)
            
            # Perform clustering to remove similar decoys
            self.log_info("Clustering similar decoys...")
            filtered_decoy_dict = self._cluster_similar_decoys(decoy_property_dict, decoy_tc_list)
            
            self.log_info(f"After clustering: {len(filtered_decoy_dict)} decoys remain")
            
            # Check if we have enough decoys
            min_decoys = self.config.param_dict['generation']['minimum_decoys_per_ligand']
            total_decoys_needed = len(lig_property_dict) * min_decoys
            
            if len(filtered_decoy_dict) < total_decoys_needed:
                self.log_error(f"Not enough decoys: need {total_decoys_needed}, have {len(filtered_decoy_dict)}")
                return False
            
            # Assign decoys to ligands using optimization
            self.log_info("Optimizing decoy assignment...")
            assignment_dict = self._assign_decoys_optimally(lig_property_dict, filtered_decoy_dict)
            
            if not assignment_dict:
                self.log_error("Decoy assignment optimization failed")
                return False
            
            # Write out final assignments
            self._write_final_assignments(lig_property_dict, filtered_decoy_dict, assignment_dict)
            
            self.log_info("Filtering step completed successfully")
            return True
            
        except Exception as e:
            self.log_error(f"Filtering step failed: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if filtering step outputs exist"""
        # Check for final assignment files
        required_files = [
            "filtered_decoys_summary.txt",
            "assignment_log.txt"
        ]
        
        for filename in required_files:
            if not os.path.exists(self.get_output_file(filename)):
                return False
        
        # Check for individual ligand assignment files
        ligand_map = self._read_ligand_map()
        for ligand_num in ligand_map.keys():
            assignment_file = self.get_output_file(f"{ligand_num}_final_property_matched_decoys.txt")
            if not os.path.exists(assignment_file):
                return False
                
        return True
    
    def _collect_ligands_and_decoys(self) -> Tuple[Dict, Dict]:
        """Collect ligands and decoys with their properties"""
        lig_property_dict = {}
        decoy_property_dict = {}
        
        # Read ligand map from setup step output
        ligand_map = self._read_ligand_map()
        
        # Collect ligand properties
        for ligand_num, (smiles, lig_id) in ligand_map.items():
            props = get_molecular_properties(smiles, get_charge=True)
            # [lig_id, smiles, mw, logp, rotb, hbd, hba, charge]
            lig_property_dict[lig_id] = [lig_id, smiles] + list(props)
        
        # Collect decoy properties
        retrieval_dir = os.path.join(os.path.dirname(self.step_dir.path), "retrieval")
        for ligand_num in ligand_map.keys():
            decoy_file = os.path.join(retrieval_dir, f"{ligand_num}_decoys.smi")
            if os.path.exists(decoy_file):
                with open(decoy_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            decoy_smiles, zinc_id = parts[0], parts[1]
                            # Get properties including charge (protonation filtering now done in retrieval)
                            props = get_molecular_properties(decoy_smiles, get_charge=True)
                            
                            # Calculate TC to closest ligand (placeholder for now)
                            tc_to_lig = 0.0
                            closest_lig = list(ligand_map.values())[0][1]  # First ligand as placeholder
                            
                            # [dec_smiles, zinc_id, mw, logp, rotb, hbd, hba, charge, prot_id, tc_to_lig, closest_lig]
                            decoy_property_dict[zinc_id] = [
                                decoy_smiles, zinc_id, props[0], props[1], props[2], 
                                props[3], props[4], props[5], "NA", tc_to_lig, closest_lig
                            ]
        
        return lig_property_dict, decoy_property_dict
    
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
    
    def _calculate_tanimoto_coefficients(self, lig_property_dict: Dict, decoy_property_dict: Dict) -> List[Tuple[float, str]]:
        """Calculate Tanimoto coefficients between decoys and ligands with caching"""
        
        # Check if we can use cached Tanimoto matrix
        tanimoto_cache_file = self.get_output_file("tanimoto_matrix.pkl")
        ligand_smiles = [props[1] for props in lig_property_dict.values()]
        decoy_smiles = [props[0] for props in decoy_property_dict.values()]
        ligand_ids = list(lig_property_dict.keys())
        decoy_ids = list(decoy_property_dict.keys())
        
        if (os.path.exists(tanimoto_cache_file) and 
            is_tanimoto_cache_valid(tanimoto_cache_file, ligand_ids, decoy_ids)):
            
            self.log_info("Loading cached Tanimoto matrix")
            try:
                tc_matrix, cached_lig_ids, cached_dec_ids = load_tanimoto_matrix(tanimoto_cache_file)
                
                # Reorder to match current order
                lig_idx_map = {lig_id: i for i, lig_id in enumerate(cached_lig_ids)}
                dec_idx_map = {dec_id: i for i, dec_id in enumerate(cached_dec_ids)}
                
            except Exception as e:
                self.log_error(f"Failed to load cached Tanimoto matrix: {e}")
                self.log_info("Recalculating Tanimoto matrix...")
                tc_matrix = calculate_tanimoto_matrix(ligand_smiles, decoy_smiles)
                save_tanimoto_matrix(tc_matrix, ligand_ids, decoy_ids, tanimoto_cache_file)
        else:
            self.log_info("Calculating Tanimoto matrix...")
            tc_matrix = calculate_tanimoto_matrix(ligand_smiles, decoy_smiles)
            save_tanimoto_matrix(tc_matrix, ligand_ids, decoy_ids, tanimoto_cache_file)
            self.log_info("Saved Tanimoto matrix to cache")
        
        # Find max TC for each decoy and update properties
        decoy_tc_list = []
        
        for j, (decoy_id, decoy_props) in enumerate(decoy_property_dict.items()):
            max_tc = 0.0
            closest_lig = None
            
            # Find maximum TC to any ligand
            for i, lig_id in enumerate(ligand_ids):
                tc = tc_matrix[i][j]
                if tc > max_tc:
                    max_tc = tc
                    closest_lig = lig_id
            
            # Update decoy properties with TC info
            decoy_property_dict[decoy_id][9] = max_tc  # tc_to_lig
            decoy_property_dict[decoy_id][10] = closest_lig  # closest_lig
            
            decoy_tc_list.append((max_tc, decoy_id))
        
        # Filter decoys within Tanimoto range
        tc_range = self.config.param_dict['property_matching']['tanimoto_range']
        min_tc, max_tc = tc_range[0], tc_range[1]
        
        filtered_tc_list = [(tc, decoy_id) for tc, decoy_id in decoy_tc_list 
                           if min_tc <= tc <= max_tc]
        
        self.log_info(f"Filtered {len(filtered_tc_list)}/{len(decoy_tc_list)} decoys within TC range [{min_tc}, {max_tc}]")
        
        return filtered_tc_list
    
    def _cluster_similar_decoys(self, decoy_property_dict: Dict, decoy_tc_list: List[Tuple[float, str]]) -> Dict:
        """Remove similar decoys using clustering (simplified version of best-first clustering)"""
        
        max_tc_between_decoys = self.config.param_dict['generation']['max_tanimoto_between_decoys']
        
        # For simplicity, we'll use a greedy approach instead of the full best-first clustering
        # This maintains the core functionality while being more manageable
        
        # Sort decoys by TC to ligands (keep most dissimilar first)
        sorted_decoys = sorted(decoy_tc_list, key=lambda x: x[0])
        
        # Generate fingerprints for all decoys
        decoy_fingerprints = {}
        for tc, decoy_id in sorted_decoys:
            if decoy_id in decoy_property_dict:
                smiles = decoy_property_dict[decoy_id][0]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    decoy_fingerprints[decoy_id] = fp
        
        # Greedy clustering: keep decoys that are sufficiently different
        cluster_heads = []
        
        for tc, decoy_id in sorted_decoys:
            if decoy_id not in decoy_fingerprints:
                continue
                
            decoy_fp = decoy_fingerprints[decoy_id]
            is_similar = False
            
            # Check similarity to existing cluster heads
            for head_id in cluster_heads:
                if head_id in decoy_fingerprints:
                    head_fp = decoy_fingerprints[head_id]
                    similarity = DataStructs.TanimotoSimilarity(decoy_fp, head_fp)
                    
                    if similarity > max_tc_between_decoys:
                        is_similar = True
                        break
            
            if not is_similar:
                cluster_heads.append(decoy_id)
        
        # Create filtered dictionary with only cluster heads
        filtered_dict = {decoy_id: decoy_property_dict[decoy_id] 
                        for decoy_id in cluster_heads 
                        if decoy_id in decoy_property_dict}
        
        self.log_info(f"Clustering reduced decoys from {len(decoy_property_dict)} to {len(filtered_dict)}")
        
        return filtered_dict
    
    def _assign_decoys_optimally(self, lig_property_dict: Dict, decoy_property_dict: Dict) -> Optional[Dict]:
        """Assign decoys to ligands using Integer Linear Programming optimization"""
        
        try:
            min_decoys = self.config.param_dict['generation']['minimum_decoys_per_ligand']
            pref_decoys = self.config.param_dict['generation']['target_decoys_per_ligand']
            
            # Create ordered lists for ILP
            lig_order = list(lig_property_dict.keys())
            decoy_order = list(decoy_property_dict.keys())
            
            # Calculate property matching windows for each decoy-ligand pair
            self.log_info("Calculating property matching windows...")
            dec_windows = []
            
            for decoy_id in decoy_order:
                decoy_props = decoy_property_dict[decoy_id]
                decoy_windows = []
                
                for lig_id in lig_order:
                    lig_props = lig_property_dict[lig_id]
                    window = self._compare_properties(lig_props, decoy_props)
                    decoy_windows.append(window)
                
                dec_windows.append(decoy_windows)
            
            # Debug: Print optimization parameters
            self.log_info(f"ILP parameters: {len(dec_windows)} decoys, {len(lig_order)} ligands, {pref_decoys} decoys per ligand")
            self.log_info(f"Total decoys needed: {len(lig_order) * pref_decoys}, available: {len(dec_windows)}")
            
            # Solve ILP assignment problem
            self.log_info("Solving optimization problem...")
            assignments = self._solve_assignment_ilp(dec_windows, pref_decoys, lig_order, decoy_order)
            
            if assignments is None:
                self.log_error("ILP optimization failed")
                return None
            
            # Convert assignments to dictionary
            assignment_dict = defaultdict(list)
            
            for i, decoy_assignment in enumerate(assignments):
                decoy_id = decoy_order[i]
                for j, assigned in enumerate(decoy_assignment):
                    if int(assigned) == 1:
                        lig_id = lig_order[j]
                        assignment_dict[lig_id].append(decoy_id)
                        break
            
            # Verify minimum requirements
            for lig_id in lig_order:
                if len(assignment_dict[lig_id]) < min_decoys:
                    self.log_error(f"Ligand {lig_id} only assigned {len(assignment_dict[lig_id])} decoys, minimum is {min_decoys}")
                    return None
            
            return dict(assignment_dict)
            
        except Exception as e:
            self.log_error(f"Assignment optimization failed: {e}")
            return None
    
    def _compare_properties(self, lig_props: List, decoy_props: List) -> Optional[int]:
        """Compare ligand and decoy properties using progressive windows"""
        
        # Extract properties as tuples
        lig_tuple = tuple(lig_props[2:8])  # (mw, logp, rotb, hbd, hba, charge)
        dec_tuple = tuple(decoy_props[2:8])
        
        # Get progressive windows
        windows = get_progressive_windows_from_config(self.config.param_dict)
        
        return compare_properties_with_windows(lig_tuple, dec_tuple, windows)
    
    def _solve_assignment_ilp(self, dec_windows: List[List[int]], pref_decoys: int, lig_order: List[str], decoy_order: List[str]) -> Optional[List[List[int]]]:
        """Solve the decoy assignment problem using Integer Linear Programming"""
        
        n = len(dec_windows)        # Number of decoys
        k = len(dec_windows[0])     # Number of ligands
        N = pref_decoys            # Target decoys per ligand
        
        # Create the problem
        prob = pulp.LpProblem("DecoyAssignmentProblem", pulp.LpMinimize)
        
        # Create decision variables only for valid assignments
        valid_assignments = []
        
        for i in range(n):
            for j in range(k):
                if dec_windows[i][j] is not None:  # Valid match
                    valid_assignments.append((i, j))
        
        if not valid_assignments:
            self.log_error("No valid decoy-ligand assignments possible")
            return None
        
        # Create decision variables only for valid assignments
        x = pulp.LpVariable.dicts("assignment", 
                                 valid_assignments,
                                 lowBound=0, upBound=1, cat=pulp.LpBinary)
        
        # Objective: minimize total matching cost (minimize "window") (only for valid assignments)
        prob += pulp.lpSum([x[(i, j)] * dec_windows[i][j] for (i, j) in valid_assignments])
        
        # Constraint: each decoy assigned to exactly one ligand
        for i in range(n):
            valid_for_decoy = [(decoy_idx, lig_idx) for (decoy_idx, lig_idx) in valid_assignments if decoy_idx == i]
            if valid_for_decoy:
                prob += pulp.lpSum([x[(decoy_idx, lig_idx)] for (decoy_idx, lig_idx) in valid_for_decoy]) == 1
        
        # Check if we can meet minimum requirements BEFORE solving
        min_decoys = self.config.param_dict['generation']['minimum_decoys_per_ligand']
        ligands_with_insufficient_decoys = []
        
        for j in range(k):
            lig_id = lig_order[j]
            this_ligand_valid = [(i, j_val) for (i, j_val) in valid_assignments if j_val == j]
            valid_count = len(this_ligand_valid)
            
            if valid_count >= min_decoys:
                # This ligand has enough valid decoys, enforce minimum
                prob += pulp.lpSum([x[(i, j_val)] for (i, j_val) in this_ligand_valid]) >= min_decoys
            else:
                # This ligand cannot meet minimum requirements
                ligands_with_insufficient_decoys.append((lig_id, valid_count, min_decoys))
        
        # If any ligands can't meet minimum, fail with helpful error message
        if ligands_with_insufficient_decoys:
            error_msg = f"Cannot meet minimum decoy requirements for {len(ligands_with_insufficient_decoys)} ligand(s):\n"
            for lig_id, available, needed in ligands_with_insufficient_decoys:
                error_msg += f"  - {lig_id}: {available} available, {needed} required\n"
            
            error_msg += "\nPossible solutions:\n"
            error_msg += f"  a) Reduce minimum_decoys_per_ligand in config (currently {min_decoys})\n"
            error_msg += f"  b) Increase max_tanimoto_between_decoys to allow more similar decoys\n"
            error_msg += f"  c) Increase total_decoys_to_generate to retrieve more options\n"
            error_msg += f"  d) Widen property_matching ranges to allow more diverse decoys\n"
            error_msg += f"  e) Check if ligands with insufficient decoys have unusual properties\n"
            
            self.log_error(error_msg)
            return None
        
        # Solve the problem
        prob.solve(pulp.PULP_CBC_CMD(msg=0))  # Suppress solver output
        
        if prob.status == 1:  # Optimal solution found
            # Convert back to the original format
            assignments = [[0 for j in range(k)] for i in range(n)]
            
            for (i, j) in valid_assignments:
                if x[(i, j)].varValue and x[(i, j)].varValue > 0.5:  # Assigned
                    assignments[i][j] = 1
            
            return assignments
        else:
            self.log_error(f"ILP solver status: {pulp.LpStatus[prob.status]}")
            return None
    
    def _apply_protonation_filtering(self, lig_property_dict: Dict, decoy_property_dict: Dict) -> Dict:
        """Apply protonation-based filtering to decoys"""
        try:
            
            # Get property windows for comparison
            windows = get_progressive_windows_from_config(self.config.param_dict)
            
            # Track which decoys match which ligands
            ligand_to_decoys = {}  # {lig_id: {decoy_id: (protomer_smiles, original_smiles, score)}}
            all_valid_decoys = {}  # Final set of decoys that match at least one ligand
            
            for lig_id, lig_props in lig_property_dict.items():
                self.log_info(f"Processing protonation filtering for ligand {lig_id}")
                
                # Extract ligand properties tuple (mw, logp, rotb, hbd, hba, charge)
                ligand_tuple = tuple(lig_props[2:8])  # Skip lig_id and smiles
                self.log_info(f"Ligand {lig_id} properties: {ligand_tuple} (charge: {ligand_tuple[5]})")
                
                # Create decoy SMILES dict for this ligand's potential decoys
                decoy_smiles_dict = {decoy_id: props[0] for decoy_id, props in decoy_property_dict.items()}
                
                try:
                    # Apply protonation filtering for this specific ligand
                    matching_decoys = self._filter_by_protonation_match(
                        ligand_tuple, 
                        decoy_smiles_dict, 
                        windows,
                        ph=7.4,
                        batch_size=50  # Process in smaller batches
                    )
                    
                    # Store ligand-specific matches
                    ligand_to_decoys[lig_id] = matching_decoys
                    
                    self.log_info(f"Ligand {lig_id}: found {len(matching_decoys)} matching decoys")
                    
                    for decoy_id, (protomer_smiles, original_smiles, score) in matching_decoys.items():
                        # Add to global valid decoys set (but preserve best protomer if multiple ligands match)
                        if decoy_id not in all_valid_decoys:
                            # Update decoy properties with protomer information
                            original_props = decoy_property_dict[decoy_id].copy()
                            
                            # Store protomer SMILES in position 8 (prot_id field)
                            original_props[8] = protomer_smiles
                            
                            # Update charge from protomer (position 7 is charge)
                            try:
                                protomer_props = get_molecular_properties(protomer_smiles, get_charge=True)
                                original_props[7] = protomer_props[5]  # Update charge from protomer
                                self.log_debug(f"Updated charge for {decoy_id}: {original_props[7]}")
                            except Exception as e:
                                self.log_warning(f"Failed to get protomer charge for {decoy_id}: {e}")
                            
                            all_valid_decoys[decoy_id] = original_props
                    
                except Exception as e:
                    self.log_error(f"Protonation filtering failed for ligand {lig_id}: {e}")
                    # Continue with other ligands
                    continue
            
            # Store ligand-to-decoys mapping for use in assignment
            self._ligand_decoy_matches = ligand_to_decoys
            
            self.log_info(f"Total unique decoys passing protonation filter: {len(all_valid_decoys)}")
            return all_valid_decoys
            
        except ImportError:
            self.log_error("Protonation module not available - skipping protonation filtering")
            self.log_error("FALLBACK: Using original decoys without protonation filtering")
            return decoy_property_dict
        except Exception as e:
            self.log_error(f"Protonation filtering failed: {e}")
            self.log_error("FALLBACK: Using original decoys without protonation filtering")
            # Fall back to original decoys if protonation fails
            return decoy_property_dict
    
    def _filter_by_protonation_match(self, ligand_props: tuple, decoy_smiles_dict: Dict[str, str], 
                                   property_windows: List[List[float]], 
                                   ph: float = 7.4, batch_size: int = 100) -> Dict[str, tuple]:
        """
        Filter decoys based on whether their protomers match ligand properties
        
        Args:
            ligand_props: (mw, logp, rotb, hbd, hba, charge) for ligand
            decoy_smiles_dict: {decoy_id: decoy_smiles} for decoys to filter
            property_windows: Progressive property tolerance windows
            ph: pH for protomer generation
            batch_size: Number of decoys to process in each batch
            
        Returns:
            Dict {decoy_id: (protomer_smiles, decoy_smiles, protomer_score)} for matching decoys
        """
        from pydock3.protonation import generate_protomers, ProtonationError
        
        if not decoy_smiles_dict:
            return {}
        
        self.log_info(f"Filtering {len(decoy_smiles_dict)} decoys based on protonation matching")
        self.log_info(f"Ligand properties: MW={ligand_props[0]:.1f}, LogP={ligand_props[1]:.2f}, "
                     f"RotB={ligand_props[2]}, HBD={ligand_props[3]}, HBA={ligand_props[4]}, Charge={ligand_props[5]}")
        
        matching_decoys = {}
        decoy_items = list(decoy_smiles_dict.items())
        
        # Process decoys in batches for efficiency
        for i in range(0, len(decoy_items), batch_size):
            batch = decoy_items[i:i + batch_size]
            self.log_debug(f"Processing batch {i//batch_size + 1}: {len(batch)} decoys")
            
            # Prepare input for protomer generation
            smiles_with_names = [f"{smiles} {decoy_id}" for decoy_id, smiles in batch]
            
            try:
                # Generate protomers for this batch
                protomer_results = generate_protomers(smiles_with_names, ph=ph)
                self.log_info(f"Generated {len(protomer_results)} total protomers for batch")
                
                # Check each protomer against ligand properties
                for protomer_smiles, decoy_id, score in protomer_results:
                    if decoy_id in decoy_smiles_dict:
                        original_smiles = decoy_smiles_dict[decoy_id]
                        
                        try:
                            # Calculate protomer properties
                            protomer_props = get_molecular_properties(protomer_smiles, get_charge=True)
                            
                            # Compare with ligand using exact charge matching
                            window = compare_properties_with_windows(ligand_props, protomer_props, property_windows)
                            
                            
                            if window is not None:
                                matching_decoys[decoy_id] = (protomer_smiles, original_smiles, score)
                                self.log_debug(f"Match found: {decoy_id} -> {protomer_smiles} (window {window}, score {score:.2f})")
                        
                        except Exception as e:
                            self.log_warning(f"Failed to process protomer {protomer_smiles} for {decoy_id}: {e}")
                            continue
                
            except ProtonationError as e:
                self.log_error(f"Batch protomer generation failed: {e}")
                # Continue with next batch rather than failing completely
                continue
        
        self.log_info(f"Found {len(matching_decoys)} decoys with matching protomers")
        return matching_decoys
    
    def _calculate_backup_decoys(self, lig_property_dict: Dict, decoy_property_dict: Dict, 
                                assignment_dict: Dict, target_decoys: int) -> Dict:
        """Split each ligand's assignments into main and backup decoys"""
        backup_dict = {}
        
        # Get progressive windows for ranking quality
        windows = get_progressive_windows_from_config(self.config.param_dict)
        
        for lig_id, assigned_decoys in assignment_dict.items():
            if len(assigned_decoys) <= target_decoys:
                # This ligand doesn't have enough for backups
                backup_dict[lig_id] = []
                continue
            
            # Get ligand properties for comparison
            lig_props = lig_property_dict[lig_id]
            lig_tuple = tuple(lig_props[2:8])  # (mw, logp, rotb, hbd, hba, charge)
            
            # Calculate quality score for each assigned decoy
            decoy_scores = []
            for decoy_id in assigned_decoys:
                decoy_props = decoy_property_dict[decoy_id]
                dec_tuple = tuple(decoy_props[2:8])
                window = compare_properties_with_windows(lig_tuple, dec_tuple, windows)
                decoy_scores.append((decoy_id, window if window is not None else 999))
            
            # Sort by quality (lower window = better match)
            decoy_scores.sort(key=lambda x: x[1])
            
            # Split into main and backup
            main_decoys = [decoy_id for decoy_id, _ in decoy_scores[:target_decoys]]
            backup_decoys = [(decoy_id, window) for decoy_id, window in decoy_scores[target_decoys:]]
            
            # Update assignment_dict to only contain main decoys
            assignment_dict[lig_id] = main_decoys
            backup_dict[lig_id] = backup_decoys
            
            self.log_info(f"Ligand {lig_id}: {len(main_decoys)} main, {len(backup_decoys)} backup decoys")
        
        return backup_dict
    
    def _write_final_assignments(self, lig_property_dict: Dict, decoy_property_dict: Dict, assignment_dict: Dict) -> None:
        """Write final decoy assignments to files with main and backup decoys"""
        
        # Get target decoys per ligand for determining main vs backup
        target_decoys = self.config.param_dict['generation']['target_decoys_per_ligand']
        
        # Calculate backup decoys for each ligand
        backup_dict = self._calculate_backup_decoys(lig_property_dict, decoy_property_dict, assignment_dict, target_decoys)
        
        # Write summary
        summary_file = self.get_output_file("filtered_decoys_summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Decoy Assignment Summary\n")
            f.write(f"========================\n")
            f.write(f"Total ligands: {len(lig_property_dict)}\n")
            f.write(f"Total decoys: {len(decoy_property_dict)}\n")
            f.write(f"Main decoys assigned: {sum(len(decoys) for decoys in assignment_dict.values())}\n")
            f.write(f"Backup decoys available: {sum(len(backups) for backups in backup_dict.values())}\n\n")
            
            for lig_id, decoys in assignment_dict.items():
                backup_count = len(backup_dict.get(lig_id, []))
                f.write(f"{lig_id}: {len(decoys)} main decoys, {backup_count} backup decoys\n")
        
        # Write assignment log (main decoys only)
        log_file = self.get_output_file("assignment_log.txt")
        with open(log_file, 'w') as f:
            f.write("Ligand_ID\tDecoy_ID\tDecoy_SMILES\tLigand_TC\tMW\tLogP\tRotB\tHBD\tHBA\tCharge\tType\n")
            
            for lig_id, decoy_ids in assignment_dict.items():
                for decoy_id in decoy_ids:
                    if decoy_id in decoy_property_dict:
                        props = decoy_property_dict[decoy_id]
                        f.write(f"{lig_id}\t{decoy_id}\t{props[0]}\t{props[9]:.3f}\t{props[2]:.1f}\t{props[3]:.2f}\t{props[4]}\t{props[5]}\t{props[6]}\t{props[7]}\tmain\n")
        
        # Write backup assignment log
        backup_log_file = self.get_output_file("backup_assignment_log.txt")
        with open(backup_log_file, 'w') as f:
            f.write("Ligand_ID\tDecoy_ID\tDecoy_SMILES\tLigand_TC\tMW\tLogP\tRotB\tHBD\tHBA\tCharge\tWindow\tRank\n")
            
            for lig_id, backup_list in backup_dict.items():
                for rank, (decoy_id, window) in enumerate(backup_list, 1):
                    if decoy_id in decoy_property_dict:
                        props = decoy_property_dict[decoy_id]
                        f.write(f"{lig_id}\t{decoy_id}\t{props[0]}\t{props[9]:.3f}\t{props[2]:.1f}\t{props[3]:.2f}\t{props[4]}\t{props[5]}\t{props[6]}\t{props[7]}\t{window}\t{rank}\n")
        
        # Write individual ligand files (matching original format)
        ligand_map = self._read_ligand_map()
        
        for ligand_num, (smiles, lig_id) in ligand_map.items():
            if lig_id in assignment_dict:
                # Get ligand properties
                lig_props = lig_property_dict[lig_id]
                
                # Main decoys file
                assignment_file = self.get_output_file(f"{ligand_num}_final_property_matched_decoys.txt")
                with open(assignment_file, 'w') as f:
                    f.write("SMILES ZINC_ID MW LogP #Rotatable_Bonds #HBond_Donors #HBond_Acceptors Charge Protomer_ID TC_TO_LIG\n")
                    
                    # Write ligand as first row
                    f.write(f"{lig_props[1]} {lig_props[0]} {lig_props[2]:.1f} {lig_props[3]:.2f} {lig_props[4]} {lig_props[5]} {lig_props[6]} {lig_props[7]} LIGAND 1.00\n")
                    
                    # Write decoys
                    for decoy_id in assignment_dict[lig_id]:
                        if decoy_id in decoy_property_dict:
                            props = decoy_property_dict[decoy_id]
                            f.write(f"{props[0]} {props[1]} {props[2]:.1f} {props[3]:.2f} {props[4]} {props[5]} {props[6]} {props[7]} {props[8]} {props[9]:.2f}\n")
                
                # Backup decoys file
                if lig_id in backup_dict and backup_dict[lig_id]:
                    backup_file = self.get_output_file(f"{ligand_num}_backup_decoys.txt")
                    with open(backup_file, 'w') as f:
                        f.write("SMILES ZINC_ID MW LogP #Rotatable_Bonds #HBond_Donors #HBond_Acceptors Charge Protomer_ID TC_TO_LIG Window Rank\n")
                        
                        # Write ligand as first row
                        f.write(f"{lig_props[1]} {lig_props[0]} {lig_props[2]:.1f} {lig_props[3]:.2f} {lig_props[4]} {lig_props[5]} {lig_props[6]} {lig_props[7]} LIGAND 1.00 0 0\n")
                        
                        # Write backup decoys
                        for rank, (decoy_id, window) in enumerate(backup_dict[lig_id], 1):
                            if decoy_id in decoy_property_dict:
                                props = decoy_property_dict[decoy_id]
                                f.write(f"{props[0]} {props[1]} {props[2]:.1f} {props[3]:.2f} {props[4]} {props[5]} {props[6]} {props[7]} {props[8]} {props[9]:.2f} {window} {rank}\n")