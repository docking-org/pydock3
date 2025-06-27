import os
import logging
from typing import Optional

from pydock3.util import Script
from pydock3.files import Dir, File
from pydock3.decoy_gen.config import DecoyGenParametersConfiguration
from pydock3.decoy_gen import __file__ as DECOY_GEN_INIT_FILE_PATH

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class DecoyGen(Script):
    JOB_DIR_NAME = "decoy_gen_job"
    CONFIG_FILE_NAME = "decoy_gen_config.yaml"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(
        os.path.dirname(DECOY_GEN_INIT_FILE_PATH),
        "default_decoy_gen_config.yaml",
    )
    WORKING_DIR_NAME = "working"

    def __init__(self):
        super().__init__()

    def new(self, job_dir_path: str = JOB_DIR_NAME) -> None:
        """Initialize new decoy generation job directory"""
        
        # Check if job dir already exists
        if os.path.exists(job_dir_path):
            logger.info(f"Job directory `{job_dir_path}` already exists. Exiting.")
            return

        # Create job dir
        job_dir = Dir(path=job_dir_path, create=True, reset=False)
        logger.info(f"Created job directory: {job_dir.path}")

        # Create working dir
        working_dir = Dir(
            path=os.path.join(job_dir.path, self.WORKING_DIR_NAME),
            create=True,
            reset=False,
        )
        logger.info(f"Created working directory: {working_dir.path}")

        # Write fresh config file from default file
        save_path = os.path.join(job_dir.path, self.CONFIG_FILE_NAME)
        DecoyGenParametersConfiguration.write_config_file(
            save_path, self.DEFAULT_CONFIG_FILE_PATH
        )
        logger.info(f"Created configuration file: {save_path}")

        # Check for actives.smi and copy it automatically
        actives_file = "actives.smi"
        if os.path.exists(actives_file):
            import shutil
            dest_path = os.path.join(job_dir.path, actives_file)
            shutil.copy2(actives_file, dest_path)
            logger.info(f"Copied {actives_file} into job directory: {dest_path}")
        
        # Check for other SMILES files in current directory
        other_smiles_files = [f for f in os.listdir(".") if f.endswith(".smi") and f != actives_file]
        if other_smiles_files:
            files_to_copy_str = "\n\t".join(other_smiles_files)
            logger.info(
                f"Found other SMILES files in current directory. Consider copying them to the job directory:\n\t{files_to_copy_str}"
            )
        elif not os.path.exists(os.path.join(job_dir.path, actives_file)):
            logger.info(
                "No SMILES files (.smi) found in current directory. You'll need to provide a SMILES file before running."
            )

    def run(
        self,
        scheduler: str,
        job_dir_path: str = ".",
        config_file_path: Optional[str] = None,
        smiles_file_path: Optional[str] = None,
        force_regenerate: bool = False,
    ) -> None:
        """Execute complete decoy generation pipeline"""
        
        logger.info(f"Starting DecoyGen pipeline with scheduler: {scheduler}")
        
        # Validate job directory
        if not os.path.exists(job_dir_path):
            raise Exception(f"Job directory does not exist: {job_dir_path}")
        
        # Determine config file path
        if config_file_path is None:
            config_file_path = os.path.join(job_dir_path, self.CONFIG_FILE_NAME)
        
        if not os.path.exists(config_file_path):
            raise Exception(f"Configuration file not found: {config_file_path}")
        
        # Load configuration
        config = DecoyGenParametersConfiguration(config_file_path)
        logger.info(f"Loaded configuration from: {config_file_path}")
        
        # Determine SMILES file path
        if smiles_file_path is None:
            smiles_file_path = os.path.join(job_dir_path, config.param_dict['input']['smiles_file'])
        
        if not os.path.exists(smiles_file_path):
            raise Exception(f"SMILES file not found: {smiles_file_path}")
        
        logger.info(f"Using SMILES file: {smiles_file_path}")
        
        # Execute pipeline steps
        working_dir_path = os.path.join(job_dir_path, self.WORKING_DIR_NAME)
        logger.info(f"Working directory: {working_dir_path}")
        logger.info(f"Building enabled: {config.param_dict['building']['enabled']}")
        
        # Run pipeline steps
        success = self._run_pipeline_steps(working_dir_path, smiles_file_path, scheduler)
        if not success:
            raise Exception("Pipeline execution failed")
        
        logger.info("DecoyGen pipeline completed successfully")
    
    def _run_pipeline_steps(self, working_dir: str, smiles_file: str, scheduler: str) -> bool:
        """Run all pipeline steps, using scheduler for retrieval if not 'local'"""
        try:
            from pydock3.decoy_gen.config import DecoyGenParametersConfiguration
            from pydock3.decoy_gen.steps.setup import SetupStep
            from pydock3.decoy_gen.steps.retrieval import RetrievalStep
            from pydock3.decoy_gen.steps.filtering import FilteringStep
            from pydock3.decoy_gen.steps.preparation import PreparationStep
            
            # Load configuration (config file is in parent directory of working dir)
            config_file = os.path.join(os.path.dirname(working_dir), self.CONFIG_FILE_NAME)
            config = DecoyGenParametersConfiguration(config_file)
            
            # Determine execution mode
            use_scheduler = scheduler != "local"
            log_func = logger.info if use_scheduler else print
            error_func = logger.error if use_scheduler else print
            
            # Run setup step (always local)
            log_func("Running Setup Step...")
            setup_step = SetupStep(working_dir, config, smiles_file)
            if not setup_step.run():
                error_func("Setup step failed")
                return False
            log_func("✓ Setup step completed")
            
            # Run retrieval step (with or without scheduler)
            retrieval_step = RetrievalStep(working_dir, config)
            if use_scheduler:
                log_func(f"Running Retrieval Step with {scheduler} scheduler...")
                retrieval_success = retrieval_step.run_with_scheduler(scheduler)
            else:
                log_func("Running Retrieval Step...")
                retrieval_success = retrieval_step.run()
                
            if not retrieval_success:
                error_func("Retrieval step failed")
                return False
            log_func("✓ Retrieval step completed")
            
            # Run filtering step (always local)
            log_func("Running Filtering Step...")
            filtering_step = FilteringStep(working_dir, config)
            if not filtering_step.run():
                error_func("Filtering step failed")
                return False
            log_func("✓ Filtering step completed")
            
            # Run preparation step (always local)
            log_func("Running Preparation Step...")
            preparation_step = PreparationStep(working_dir, config)
            if not preparation_step.run():
                error_func("Preparation step failed")
                return False
            log_func("✓ Preparation step completed")
            
            # TODO: Implement building steps if enabled
            if config.param_dict['building']['enabled']:
                log_func("Building steps are enabled but not yet implemented")
            
            return True
            
        except Exception as e:
            error_func = logger.error if scheduler != "local" else print
            error_func(f"Pipeline execution failed: {e}")
            import traceback
            traceback.print_exc()
            return False