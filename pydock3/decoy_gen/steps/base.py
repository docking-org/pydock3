import os
import logging
from abc import ABC, abstractmethod
from typing import Optional

from pydock3.decoy_gen.config import DecoyGenParametersConfiguration
from pydock3.files import Dir, File

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class DecoyGenStep(ABC):
    """Base class for decoy generation steps"""
    
    def __init__(self, working_dir: str, config: DecoyGenParametersConfiguration, step_name: str):
        self.working_dir = working_dir
        self.config = config
        self.step_name = step_name
        self.step_dir = Dir(
            path=os.path.join(working_dir, step_name),
            create=True,
            reset=False
        )
        
    @abstractmethod
    def run(self) -> bool:
        """Execute the step"""
        pass
        
    @abstractmethod
    def is_complete(self) -> bool:
        """Check if step outputs exist"""
        pass
        
    def cleanup(self) -> None:
        """Clean intermediate files - default implementation does nothing"""
        logger.debug(f"Cleanup called for {self.step_name} - no cleanup implemented")
        
    def get_output_file(self, filename: str) -> str:
        """Get path to output file in step directory"""
        return os.path.join(self.step_dir.path, filename)
        
    def log_info(self, message: str) -> None:
        """Log info message with step name"""
        logger.info(f"[{self.step_name}] {message}")
        
    def log_debug(self, message: str) -> None:
        """Log debug message with step name"""
        logger.debug(f"[{self.step_name}] {message}")
        
    def log_error(self, message: str) -> None:
        """Log error message with step name"""
        logger.error(f"[{self.step_name}] {message}")


class DecoyGenPipeline:
    """Manages the sequence of decoy generation steps"""
    
    def __init__(self, working_dir: str, config: DecoyGenParametersConfiguration):
        self.working_dir = working_dir
        self.config = config
        self.steps = []
        
    def add_step(self, step: DecoyGenStep) -> None:
        """Add a step to the pipeline"""
        self.steps.append(step)
        
    def run_all(self, force_regenerate: bool = False) -> bool:
        """Run all steps in sequence"""
        logger.info(f"Starting decoy generation pipeline with {len(self.steps)} steps")
        
        for i, step in enumerate(self.steps, 1):
            logger.info(f"Step {i}/{len(self.steps)}: {step.step_name}")
            
            # Check if step is already complete
            if not force_regenerate and step.is_complete():
                logger.info(f"Step {step.step_name} already complete, skipping")
                continue
                
            # Run the step
            success = step.run()
            if not success:
                logger.error(f"Step {step.step_name} failed")
                return False
                
            # Verify completion
            if not step.is_complete():
                logger.error(f"Step {step.step_name} did not produce expected outputs")
                return False
                
        logger.info("All pipeline steps completed successfully")
        return True