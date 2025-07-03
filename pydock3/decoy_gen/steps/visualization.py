import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import numpy as np
import logging

# Suppress matplotlib debug messages while keeping other debug logs
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

from pydock3.decoy_gen.steps.base import DecoyGenStep

class VisualizationStep(DecoyGenStep):
    """Visualization step: creates plots comparing ligand and decoy property distributions"""
    
    def __init__(self, working_dir: str, config):
        super().__init__(working_dir, config, "visualization")
        
    def run(self) -> bool:
        """Execute visualization step"""
        try:
            self.log_info("Starting visualization step")
            
            # Verify prerequisite step outputs exist
            preparation_dir = os.path.join(os.path.dirname(self.step_dir.path), "preparation")
            all_decoys_file = os.path.join(preparation_dir, "all_decoys.smi")
            
            if not os.path.exists(all_decoys_file):
                self.log_error("Preparation step must be completed before visualization")
                return False
            
            # Load data
            ligand_data, decoy_data, assignments = self._load_data()
            
            if ligand_data.empty or decoy_data.empty:
                self.log_error("No data found for visualization")
                return False
            
            self.log_info(f"Loaded {len(ligand_data)} ligands and {len(decoy_data)} decoys")
            
            # Create overall property comparison plots
            self._create_overall_comparison_plots(ligand_data, decoy_data)
            
            # Create per-ligand plots
            self._create_per_ligand_plots(ligand_data, decoy_data, assignments)
            
            self.log_info("Visualization step completed successfully")
            return True
            
        except Exception as e:
            self.log_error(f"Visualization step failed: {str(e)}")
            import traceback
            traceback.print_exc()
            return False
    
    def is_complete(self) -> bool:
        """Check if visualization step outputs exist"""
        required_files = [
            "overall_property_comparison.png",
            "property_distribution_summary.txt"
        ]
        
        for filename in required_files:
            if not os.path.exists(self.get_output_file(filename)):
                return False
                
        return True
    
    def _load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
        """Load ligand and decoy data from previous steps"""
        
        # Load ligand data from setup step
        setup_dir = os.path.join(os.path.dirname(self.step_dir.path), "setup")
        ligand_map_file = os.path.join(setup_dir, "LIGAND_MAP.txt")
        
        ligand_data = []
        with open(ligand_map_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    ligand_num, smiles, lig_id = parts[0], parts[1], parts[2]
                    ligand_data.append({
                        'ligand_num': ligand_num,
                        'smiles': smiles,
                        'ligand_id': lig_id
                    })
        
        # Load decoy data from preparation step
        preparation_dir = os.path.join(os.path.dirname(self.step_dir.path), "preparation")
        all_decoys_file = os.path.join(preparation_dir, "all_decoys.smi")
        
        decoy_data = []
        assignments = {}
        
        with open(all_decoys_file, 'r') as f:
            header = f.readline().strip().split()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 10:
                    decoy_info = {
                        'smiles': parts[0],
                        'zinc_id': parts[1],
                        'ligand_id': parts[2],
                        'tc_to_lig': float(parts[3]),
                        'mw': float(parts[4]),
                        'logp': float(parts[5]),
                        'rotb': int(parts[6]),
                        'hbd': int(parts[7]),
                        'hba': int(parts[8]),
                        'charge': int(parts[9])
                    }
                    decoy_data.append(decoy_info)
                    
                    # Build assignments dictionary
                    lig_id = decoy_info['ligand_id']
                    if lig_id not in assignments:
                        assignments[lig_id] = []
                    assignments[lig_id].append(decoy_info)
        
        # Calculate ligand properties and add to ligand_data
        from pydock3.decoy_gen.utils import get_molecular_properties
        
        for lig_info in ligand_data:
            props = get_molecular_properties(lig_info['smiles'], get_charge=True)
            lig_info.update({
                'mw': props[0],
                'logp': props[1], 
                'rotb': props[2],
                'hbd': props[3],
                'hba': props[4],
                'charge': props[5]
            })
        
        return pd.DataFrame(ligand_data), pd.DataFrame(decoy_data), assignments
    
    def _create_overall_comparison_plots(self, ligand_data: pd.DataFrame, decoy_data: pd.DataFrame):
        """Create overall comparison plots for all properties"""
        
        properties = [
            ('mw', 'Molecular Weight (Da)', 'continuous'),
            ('logp', 'LogP', 'continuous'),
            ('rotb', 'Rotatable Bonds', 'discrete'),
            ('hbd', 'H-Bond Donors', 'discrete'),
            ('hba', 'H-Bond Acceptors', 'discrete'),
            ('charge', 'Formal Charge', 'discrete')
        ]
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('Ligand vs Decoy Property Distributions', fontsize=16, fontweight='bold')
        
        axes = axes.flatten()
        
        for idx, (prop, label, plot_type) in enumerate(properties):
            ax = axes[idx]
            
            if plot_type == 'continuous':
                # Use step plots instead of histograms for better comparison with different sample sizes
                lig_vals = ligand_data[prop].values
                dec_vals = decoy_data[prop].values
                
                # Create bins covering the full range
                all_vals = np.concatenate([lig_vals, dec_vals])
                bins = np.linspace(all_vals.min(), all_vals.max(), 30)
                
                # Calculate normalized histograms (percentage of total)
                lig_counts, _ = np.histogram(lig_vals, bins=bins)
                dec_counts, _ = np.histogram(dec_vals, bins=bins)
                
                lig_pct = lig_counts / len(lig_vals) * 100
                dec_pct = dec_counts / len(dec_vals) * 100
                
                # Plot as step functions for cleaner comparison
                bin_centers = (bins[:-1] + bins[1:]) / 2
                ax.step(bin_centers, lig_pct, where='mid', linewidth=2, alpha=0.8, 
                       label=f'Ligands (n={len(ligand_data)})', color='blue')
                ax.step(bin_centers, dec_pct, where='mid', linewidth=2, alpha=0.7, 
                       label=f'Decoys (n={len(decoy_data)})', color='orange')
                
                ax.set_ylabel('Percentage (%)')
            else:
                # Bar plot for discrete properties
                lig_counts = ligand_data[prop].value_counts().sort_index()
                dec_counts = decoy_data[prop].value_counts().sort_index()
                
                # Normalize to percentages
                lig_pct = lig_counts / len(ligand_data) * 100
                dec_pct = dec_counts / len(decoy_data) * 100
                
                x_vals = sorted(set(ligand_data[prop].tolist() + decoy_data[prop].tolist()))
                
                lig_vals = [lig_pct.get(x, 0) for x in x_vals]
                dec_vals = [dec_pct.get(x, 0) for x in x_vals]
                
                x_pos = np.arange(len(x_vals))
                width = 0.35
                
                ax.bar(x_pos - width/2, lig_vals, width, alpha=0.7, label=f'Ligands (n={len(ligand_data)})', color='blue')
                ax.bar(x_pos + width/2, dec_vals, width, alpha=0.5, label=f'Decoys (n={len(decoy_data)})', color='orange')
                ax.set_xticks(x_pos)
                ax.set_xticklabels(x_vals)
                ax.set_ylabel('Percentage (%)')
            
            ax.set_xlabel(label)
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Add statistics
            lig_mean = ligand_data[prop].mean()
            dec_mean = decoy_data[prop].mean()
            ax.text(0.02, 0.98, f'Ligand μ={lig_mean:.2f}\nDecoy μ={dec_mean:.2f}', 
                   transform=ax.transAxes, verticalalignment='top', 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(self.get_output_file("overall_property_comparison.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        self.log_info("Created overall property comparison plot")
    
    def _create_per_ligand_plots(self, ligand_data: pd.DataFrame, decoy_data: pd.DataFrame, assignments: Dict):
        """Create individual plots for each ligand showing assigned decoy distributions"""
        
        properties = [
            ('mw', 'Molecular Weight (Da)'),
            ('logp', 'LogP'),
            ('rotb', 'Rotatable Bonds'),
            ('hbd', 'H-Bond Donors'),
            ('hba', 'H-Bond Acceptors'),
            ('charge', 'Formal Charge')
        ]
        
        summary_stats = []
        
        for _, ligand in ligand_data.iterrows():
            lig_id = ligand['ligand_id']
            
            if lig_id not in assignments or len(assignments[lig_id]) == 0:
                continue
            
            # Get assigned decoys for this ligand
            assigned_decoys = pd.DataFrame(assignments[lig_id])
            
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f'Ligand {lig_id} - Assigned Decoy Distributions\nLigand SMILES: {ligand["smiles"]}', 
                        fontsize=14, fontweight='bold')
            
            axes = axes.flatten()
            
            for idx, (prop, label) in enumerate(properties):
                ax = axes[idx]
                
                # Plot histogram of assigned decoys
                decoy_vals = assigned_decoys[prop]
                lig_val = ligand[prop]
                
                if len(decoy_vals.unique()) > 5:  # Continuous-like
                    ax.hist(decoy_vals, bins=min(15, len(decoy_vals)//2 + 1), alpha=0.7, 
                           color='lightblue', edgecolor='black', label=f'Assigned Decoys (n={len(decoy_vals)})')
                    ax.axvline(lig_val, color='red', linewidth=3, label=f'Ligand ({lig_val})', linestyle='--')
                else:  # Discrete
                    counts = decoy_vals.value_counts().sort_index()
                    ax.bar(counts.index, counts.values, alpha=0.7, color='lightblue', edgecolor='black')
                    ax.axvline(lig_val, color='red', linewidth=3, label=f'Ligand ({lig_val})', linestyle='--')
                
                ax.set_xlabel(label)
                ax.set_ylabel('Count')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                # Add statistics
                decoy_mean = decoy_vals.mean()
                decoy_std = decoy_vals.std()
                ax.text(0.02, 0.98, f'Decoy μ={decoy_mean:.2f}\nDecoy σ={decoy_std:.2f}', 
                       transform=ax.transAxes, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Store summary statistics
                summary_stats.append({
                    'ligand_id': lig_id,
                    'property': prop,
                    'ligand_value': lig_val,
                    'decoy_mean': decoy_mean,
                    'decoy_std': decoy_std,
                    'decoy_min': decoy_vals.min(),
                    'decoy_max': decoy_vals.max(),
                    'n_decoys': len(decoy_vals)
                })
            
            plt.tight_layout()
            plt.savefig(self.get_output_file(f"ligand_{lig_id}_property_distributions.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            self.log_info(f"Created property distribution plot for ligand {lig_id}")
        
        # Write summary statistics
        self._write_distribution_summary(summary_stats)
    
    def _write_distribution_summary(self, summary_stats: List[Dict]):
        """Write summary statistics to file"""
        
        summary_file = self.get_output_file("property_distribution_summary.txt")
        
        with open(summary_file, 'w') as f:
            f.write("Property Distribution Summary\n")
            f.write("=" * 50 + "\n\n")
            
            # Group by ligand
            ligands = {}
            for stat in summary_stats:
                lig_id = stat['ligand_id']
                if lig_id not in ligands:
                    ligands[lig_id] = []
                ligands[lig_id].append(stat)
            
            for lig_id, stats in ligands.items():
                f.write(f"Ligand: {lig_id}\n")
                f.write("-" * 30 + "\n")
                f.write(f"{'Property':<15} {'Ligand':<8} {'Decoy Mean':<12} {'Decoy Std':<10} {'Range':<15} {'N':<5}\n")
                f.write("-" * 70 + "\n")
                
                for stat in stats:
                    prop = stat['property'].upper()
                    lig_val = stat['ligand_value']
                    dec_mean = stat['decoy_mean']
                    dec_std = stat['decoy_std']
                    dec_range = f"{stat['decoy_min']:.1f}-{stat['decoy_max']:.1f}"
                    n_decoys = stat['n_decoys']
                    
                    f.write(f"{prop:<15} {lig_val:<8.2f} {dec_mean:<12.2f} {dec_std:<10.2f} {dec_range:<15} {n_decoys:<5}\n")
                
                f.write("\n")
        
        self.log_info("Created property distribution summary")