import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass

from ucsfdock.util import system_call

from ucsfdock.submit.slurm import __file__ as DOCK_SUBMISSION_SLURM_INIT_FILE_PATH
DOCK_SUBMISSION_SLURM_DIR_PATH = os.path.dirname(DOCK_SUBMISSION_SLURM_INIT_FILE_PATH)

from ucsfdock.submit.sge import __file__ as DOCK_SUBMISSION_SGE_INIT_FILE_PATH
DOCK_SUBMISSION_SGE_DIR_PATH = os.path.dirname(DOCK_SUBMISSION_SGE_INIT_FILE_PATH)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@dataclass
class JobScheduler(ABC):
    name: str

    def __post_init__(self):
        #
        if (os.environ.get("SHRTCACHE") is None) or (os.environ.get("LONGCACHE") is None):
            raise Exception(
                "Need to set environmental variables SHRTCACHE and LONGCACHE in order to use a job scheduler.")
        self.SHRTCACHE = os.environ["SHRTCACHE"]
        self.LONGCACHE = os.environ["LONGCACHE"]

    @abstractmethod
    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        raise NotImplementedError

    @abstractmethod
    def is_running_job(self, job_name):
        raise NotImplementedError


@dataclass
class NoJobScheduler(JobScheduler):
    name: str = "None"

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        raise NotImplementedError

    def is_running_job(self, job_name):
        raise NotImplementedError


@dataclass
class SlurmJobScheduler(JobScheduler):
    name: str = "Slurm"
    submissision_script_path: str = os.path.join(DOCK_SUBMISSION_SLURM_DIR_PATH, "subdock.bash")

    def __post_init__(self):
        #
        self.SBATCH_EXEC = os.environ.get("SBATCH_EXEC")
        self.SLURM_SETTINGS = os.environ.get("SLURM_SETTINGS")
        self.SQUEUE_EXEC = os.environ.get("SQUEUE_EXEC")
        required_env_var_names = [
            "SQUEUE_EXEC",
            "SBATCH_EXEC",
        ]
        if any([globals()[env_var] is None for env_var in required_env_var_names]):
            raise Exception(f"The following environmental variables are required to use the Slurm job scheduler: {required_env_var_names}")

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        command_str = f"bash {self.submissision_script_path}"
        env_vars_dict = {
            "EXPORT_DEST": output_dir.path,
            "INPUT_SOURCE": input_file.path,
            "DOCKFILES": dock_files_dir.path,
            "DOCKEXEC": dock_executable_path,
            "SHRTCACHE": self.SHRTCACHE,
            "LONGCACHE": self.LONGCACHE,
        }

        if self.SBATCH_EXEC is not None:
            env_vars_dict["SBATCH_EXEC"] = self.SBATCH_EXEC
        if self.SLURM_SETTINGS is not None:
            env_vars_dict["SLURM_SETTINGS"] = self.SLURM_SETTINGS

        sbatch_args_str = f"-J {job_name} -o {output_dir.path}/{job_name}_%A_%a.out -e {output_dir.path}/{job_name}_%A_%a.err"

        if job_timeout_minutes is None:
            sbatch_args_str += "--time=0 "
        else:
            sbatch_args_str += f"--time={job_timeout_minutes} "

        env_vars_dict["SBATCH_ARGS"] = sbatch_args_str

        return system_call(command_str, cwd=working_dir.path, env_vars_dict=env_vars_dict)       

    def is_running_job(self, job_name):
        command_str = f"{self.SQUEUE_EXEC} --format='%.18i %.{len(job_name)}j' | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False


@dataclass
class SGEJobScheduler(JobScheduler):
    name: str = "SGE"
    submissision_script_path: str = os.path.join(DOCK_SUBMISSION_SGE_DIR_PATH, "subdock.bash")

    def __post_init__(self):
        #
        self.QSUB_EXEC = os.environ.get("QSUB_EXEC")
        self.SGE_SETTINGS = os.environ.get("SGE_SETTINGS")
        self.QSTAT_EXEC = os.environ.get("QSTAT_EXEC")
        required_env_var_names = [
            "QSTAT_EXEC",
            "QSUB_EXEC",
        ]
        if any([globals()[env_var] is None for env_var in required_env_var_names]):
            raise Exception(f"The following environmental variables are required to use the SGE job scheduler: {required_env_var_names}")

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        command_str = f"bash {self.submissision_script_path}"
        env_vars_dict = {
            "EXPORT_DEST": output_dir.path,
            "INPUT_SOURCE": input_file.path,
            "DOCKFILES": dock_files_dir.path,
            "DOCKEXEC": dock_executable_path,
            "SHRTCACHE": self.SHRTCACHE,
            "LONGCACHE": self.LONGCACHE,
        }

        if self.QSUB_EXEC is not None:
            env_vars_dict["QSUB_EXEC"] = self.QSUB_EXEC
        if self.SGE_SETTINGS is not None:
            env_vars_dict["SGE_SETTINGS"] = self.SGE_SETTINGS

        qsub_args_str = f"-N {job_name} -o {output_dir.path} -e {output_dir.path}"

        if job_timeout_minutes is not None:
            job_timeout_seconds = 60 * job_timeout_minutes
            qsub_args_str += f"-l s_rt={job_timeout_seconds} -l h_rt={job_timeout_seconds} "

        env_vars_dict["QSUB_ARGS"] = qsub_args_str

        return system_call(command_str, cwd=working_dir.path, env_vars_dict=env_vars_dict)

    def is_running_job(self, job_name):
        command_str = f"{self.QSTAT_EXEC} -r | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False
