"""
Job manager for decoy generation retrieval step
"""

import os
import time
import logging
import shutil
from typing import List, Dict, Set, Tuple

from pydock3.job_schedulers import SlurmJobScheduler, SGEJobScheduler
from pydock3.decoy_gen.steps.retrieval import SingleLigandRetrievalStep


logger = logging.getLogger(__name__)


class RetrievalJobManager:
    """Manages job submission and monitoring for ligand retrieval using submit_single_step"""
    
    def __init__(self, working_dir: str, config, scheduler_name: str):
        self.working_dir = working_dir
        self.config = config
        self.scheduler = self._get_scheduler(scheduler_name)
        
        # Job configuration from config
        job_config = config.param_dict.get('job_submission', {})
        self.max_retries = job_config.get('max_retries_per_ligand', 3)
        self.job_timeout = job_config.get('job_timeout_minutes', 60)
        
        # Job tracking
        self.ligand_retry_count: Dict[str, int] = {}
        self.completed_ligands: Set[str] = set()
        self.failed_ligands: Set[str] = set()
        self.submitted_jobs: Dict[str, str] = {}  # ligand_num -> job_name
        
    def _get_scheduler(self, scheduler_name: str):
        """Get appropriate scheduler instance"""
        if scheduler_name.lower() == 'slurm':
            return SlurmJobScheduler()
        elif scheduler_name.lower() == 'sge':
            return SGEJobScheduler()
        else:
            raise ValueError(f"Unsupported scheduler: {scheduler_name}")
    
    def submit_and_wait_for_ligands(self, ligand_map: Dict[str, Tuple[str, str]]) -> bool:
        """Submit jobs for all ligands and wait for completion"""
        
        # Initialize retry counts
        for ligand_num in ligand_map.keys():
            self.ligand_retry_count[ligand_num] = 0
        
        # Submit initial jobs for all ligands
        pending_ligands = set(ligand_map.keys())
        self._submit_jobs_for_ligands(ligand_map, pending_ligands)
        
        # Monitor and retry failed jobs
        while pending_ligands:
            time.sleep(30)  # Check every 30 seconds
            
            newly_completed, newly_failed = self._check_job_status(pending_ligands)
            
            # Handle completed jobs
            for ligand_num in newly_completed:
                self.completed_ligands.add(ligand_num)
                pending_ligands.discard(ligand_num)
                logger.info(f"✓ Ligand {ligand_num} retrieval completed")
                
                # Copy result to main retrieval directory
                self._copy_result_to_main_dir(ligand_num)
            
            # Handle failed jobs
            retry_ligands = {}
            for ligand_num in newly_failed:
                self.ligand_retry_count[ligand_num] += 1
                
                if self.ligand_retry_count[ligand_num] <= self.max_retries:
                    retry_ligands[ligand_num] = ligand_map[ligand_num]
                    logger.warning(f"⚠ Ligand {ligand_num} failed, retrying ({self.ligand_retry_count[ligand_num]}/{self.max_retries})")
                else:
                    self.failed_ligands.add(ligand_num)
                    pending_ligands.discard(ligand_num)
                    logger.error(f"✗ Ligand {ligand_num} failed after {self.max_retries} retries, skipping")
            
            # Resubmit failed jobs
            if retry_ligands:
                self._submit_jobs_for_ligands(ligand_map, retry_ligands.keys())
        
        # Report results
        total_ligands = len(ligand_map)
        completed_count = len(self.completed_ligands)
        failed_count = len(self.failed_ligands)
        
        logger.info(f"Retrieval completed: {completed_count}/{total_ligands} ligands successful, {failed_count} failed")
        
        # Return success if at least some ligands completed
        return completed_count > 0
    
    def _submit_jobs_for_ligands(self, ligand_map: Dict[str, Tuple[str, str]], ligand_nums: List[str]) -> None:
        """Submit jobs for specified ligands using submit_single_step"""
        
        for ligand_num in ligand_nums:
            smiles, lig_id = ligand_map[ligand_num]
            
            # Create step instance for this ligand
            step_instance = SingleLigandRetrievalStep(
                working_dir=self.working_dir,
                config=self.config,
                ligand_num=ligand_num,
                ligand_smiles=smiles,
                ligand_id=lig_id
            )
            
            # Generate unique job name
            job_name = f"decoy_retrieval_{ligand_num}_{int(time.time())}"
            self.submitted_jobs[ligand_num] = job_name
            
            # Submit using existing scheduler method
            try:
                self.scheduler.submit_single_step(
                    step_instance=step_instance,
                    job_name=job_name
                )
                logger.info(f"Submitted job {job_name} for ligand {ligand_num}")
            except Exception as e:
                logger.error(f"Failed to submit job for ligand {ligand_num}: {e}")
                # Mark as failed immediately
                self.failed_ligands.add(ligand_num)
        
        logger.info(f"Submitted jobs for ligands: {list(ligand_nums)}")
    
    def _check_job_status(self, pending_ligands: Set[str]) -> Tuple[List[str], List[str]]:
        """Check status of pending jobs using scheduler's job_is_on_queue"""
        completed = []
        failed = []
        
        for ligand_num in list(pending_ligands):
            job_name = self.submitted_jobs.get(ligand_num)
            if not job_name:
                continue
                
            # Check if job is still running
            try:
                if self.scheduler.job_is_on_queue(job_name):
                    continue  # Still running
            except Exception as e:
                logger.warning(f"Error checking job status for {job_name}: {e}")
                continue
            
            # Job finished, check if it completed successfully
            step_dir = os.path.join(self.working_dir, "retrieval", f"job_{ligand_num}")
            output_file = os.path.join(step_dir, f"{ligand_num}_decoys.smi")
            
            if os.path.exists(output_file):
                completed.append(ligand_num)
            else:
                failed.append(ligand_num)
        
        return completed, failed
    
    def _copy_result_to_main_dir(self, ligand_num: str) -> None:
        """Copy result from job directory to main retrieval directory"""
        source_dir = os.path.join(self.working_dir, "retrieval", f"job_{ligand_num}")
        dest_dir = os.path.join(self.working_dir, "retrieval")
        
        source_file = os.path.join(source_dir, f"{ligand_num}_decoys.smi")
        dest_file = os.path.join(dest_dir, f"{ligand_num}_decoys.smi")
        
        if os.path.exists(source_file):
            os.makedirs(dest_dir, exist_ok=True)
            try:
                shutil.copy2(source_file, dest_file)
                logger.info(f"Copied result for ligand {ligand_num} to main directory")
            except Exception as e:
                logger.error(f"Failed to copy result for ligand {ligand_num}: {e}")
        else:
            logger.error(f"Expected result file not found: {source_file}")