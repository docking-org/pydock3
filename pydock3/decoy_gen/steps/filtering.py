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
            
            # Calculate Tanimoto coefficients if enabled
            if self.config.param_dict['generation']['calculate_tanimoto']:
                self.log_info("Calculating Tanimoto coefficients...")
                decoy_tc_list = self._calculate_tanimoto_coefficients(lig_property_dict, decoy_property_dict)
            else:
                self.log_info("Tanimoto calculation disabled")
                decoy_tc_list = [(0.0, decoy_id) for decoy_id in decoy_property_dict.keys()]
            
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
            props = get_molecular_properties(smiles)
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
                            props = get_molecular_properties(decoy_smiles)
                            
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
                    decoy_windows.append(window if window is not None else 1000)  # 1000 = bad match
                
                dec_windows.append(decoy_windows)
            
            # Debug: Print optimization parameters
            self.log_info(f"ILP parameters: {len(dec_windows)} decoys, {len(lig_order)} ligands, {pref_decoys} decoys per ligand")
            self.log_info(f"Total decoys needed: {len(lig_order) * pref_decoys}, available: {len(dec_windows)}")
            
            # Solve ILP assignment problem
            self.log_info("Solving optimization problem...")
            assignments = self._solve_assignment_ilp(dec_windows, pref_decoys)
            
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
    
    def _solve_assignment_ilp(self, dec_windows: List[List[int]], pref_decoys: int) -> Optional[List[List[int]]]:
        """Solve the decoy assignment problem using Integer Linear Programming"""
        
        n = len(dec_windows)        # Number of decoys
        k = len(dec_windows[0])     # Number of ligands
        N = pref_decoys            # Target decoys per ligand
        
        # Create the problem
        prob = pulp.LpProblem("DecoyAssignmentProblem", pulp.LpMinimize)
        
        # Create decision variables
        x = pulp.LpVariable.dicts("assignment", 
                                 ((i, j) for i in range(n) for j in range(k)),
                                 lowBound=0, upBound=1, cat=pulp.LpBinary)
        
        # Objective: minimize total matching cost
        prob += pulp.lpSum([x[(i, j)] * dec_windows[i][j] for i in range(n) for j in range(k)])
        
        # Constraint: each decoy assigned to exactly one ligand
        for i in range(n):
            prob += pulp.lpSum([x[(i, j)] for j in range(k)]) == 1
        
        # Constraint: each ligand gets at least minimum decoys (use minimum, not preferred)
        min_decoys = self.config.param_dict['generation']['minimum_decoys_per_ligand']
        for j in range(k):
            prob += pulp.lpSum([x[(i, j)] for i in range(n)]) >= min_decoys
        
        # Solve the problem
        prob.solve(pulp.PULP_CBC_CMD(msg=0))  # Suppress solver output
        
        if prob.status == 1:  # Optimal solution found
            assignments = [[x[(i, j)].varValue for j in range(k)] for i in range(n)]
            return assignments
        else:
            self.log_error(f"ILP solver status: {pulp.LpStatus[prob.status]}")
            return None
    
    def _write_final_assignments(self, lig_property_dict: Dict, decoy_property_dict: Dict, assignment_dict: Dict) -> None:
        """Write final decoy assignments to files"""
        
        # Write summary
        summary_file = self.get_output_file("filtered_decoys_summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Decoy Assignment Summary\n")
            f.write(f"========================\n")
            f.write(f"Total ligands: {len(lig_property_dict)}\n")
            f.write(f"Total decoys: {len(decoy_property_dict)}\n")
            f.write(f"Decoys assigned: {sum(len(decoys) for decoys in assignment_dict.values())}\n\n")
            
            for lig_id, decoys in assignment_dict.items():
                f.write(f"{lig_id}: {len(decoys)} decoys assigned\n")
        
        # Write assignment log
        log_file = self.get_output_file("assignment_log.txt")
        with open(log_file, 'w') as f:
            f.write("Ligand_ID\tDecoy_ID\tDecoy_SMILES\tLigand_TC\tMW\tLogP\tRotB\tHBD\tHBA\tCharge\n")
            
            for lig_id, decoy_ids in assignment_dict.items():
                for decoy_id in decoy_ids:
                    if decoy_id in decoy_property_dict:
                        props = decoy_property_dict[decoy_id]
                        f.write(f"{lig_id}\t{decoy_id}\t{props[0]}\t{props[9]:.3f}\t{props[2]:.1f}\t{props[3]:.2f}\t{props[4]}\t{props[5]}\t{props[6]}\t{props[7]}\n")
        
        # Write individual ligand files (matching original format)
        ligand_map = self._read_ligand_map()
        
        for ligand_num, (smiles, lig_id) in ligand_map.items():
            if lig_id in assignment_dict:
                assignment_file = self.get_output_file(f"{ligand_num}_final_property_matched_decoys.txt")
                with open(assignment_file, 'w') as f:
                    f.write("SMILES ZINC_ID logP #Rotatable_Bonds #HBond_Donors #HBond_Acceptors Charge Protomer_ID TC_TO_LIG\n")
                    
                    for decoy_id in assignment_dict[lig_id]:
                        if decoy_id in decoy_property_dict:
                            props = decoy_property_dict[decoy_id]
                            f.write(f"{props[0]} {props[1]} {props[3]} {props[4]} {props[5]} {props[6]} {props[7]} {props[8]} {props[9]:.6f}\n")