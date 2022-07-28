import logging
import os
from abc import ABC, abstractmethod

from ucsfdock.util import system_call
from ucsfdock.submit.slurm import __file__ as DOCK_SUBMISSION_SLURM_INIT_FILE_PATH
from ucsfdock.submit.sge import __file__ as DOCK_SUBMISSION_SGE_INIT_FILE_PATH

DOCK_SUBMISSION_SLURM_DIR_PATH = os.path.dirname(DOCK_SUBMISSION_SLURM_INIT_FILE_PATH)
DOCK_SUBMISSION_SGE_DIR_PATH = os.path.dirname(DOCK_SUBMISSION_SGE_INIT_FILE_PATH)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class JobScheduler(ABC):
    SUBMISSION_SCRIPT_PATH = None
    REQUIRED_ENV_VAR_NAMES = []

    def __init__(self, name):
        self.name = name

    @abstractmethod
    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        raise NotImplementedError

    @abstractmethod
    def is_running_job(self, job_name):
        raise NotImplementedError


class NoJobScheduler(JobScheduler):
    SUBMISSION_SCRIPT_PATH = None
    REQUIRED_ENV_VAR_NAMES = []

    def __init__(self, name):
        super().__init__(name)

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        raise NotImplementedError

    def is_running_job(self, job_name):
        raise NotImplementedError


class SlurmJobScheduler(JobScheduler):
    SUBMISSION_SCRIPT_PATH = os.path.join(DOCK_SUBMISSION_SLURM_DIR_PATH, "subdock.bash")
    REQUIRED_ENV_VAR_NAMES = [
        "SHRTCACHE",
        "LONGCACHE",
        "SBATCH_EXEC",
        "SQUEUE_EXEC",
    ]

    def __init__(self):
        super().__init__(name="Slurm")

        # set required env vars
        self.SHRTCACHE = os.environ["SHRTCACHE"]
        self.LONGCACHE = os.environ["LONGCACHE"]
        self.SBATCH_EXEC = os.environ["SBATCH_EXEC"]
        self.SQUEUE_EXEC = os.environ["SQUEUE_EXEC"]

        # set optional env vars
        self.SLURM_SETTINGS = os.environ.get("SLURM_SETTINGS")

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        command_str = f"bash {self.SUBMISSION_SCRIPT_PATH}"
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

        sbatch_args_str = f"-J {job_name} -o {output_dir.path}/{job_name}_%A_%a.out -e {output_dir.path}/{job_name}_%A_%a.err "

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
        

class SGEJobScheduler(JobScheduler):
    SUBMISSION_SCRIPT_PATH = os.path.join(DOCK_SUBMISSION_SGE_DIR_PATH, "subdock.bash")
    REQUIRED_ENV_VAR_NAMES = [
        "SHRTCACHE",
        "LONGCACHE",
        "QSUB_EXEC",
        "QSTAT_EXEC",
    ]

    def __init__(self):
        super().__init__(name="SGE")

        # set required env vars
        self.SHRTCACHE = os.environ["SHRTCACHE"]
        self.LONGCACHE = os.environ["LONGCACHE"]
        self.QSUB_EXEC = os.environ["QSUB_EXEC"]
        self.QSTAT_EXEC = os.environ["QSTAT_EXEC"]

        # set optional env vars
        self.SGE_SETTINGS = os.environ.get("SGE_SETTINGS")

    def run(self, job_name, dock_executable_path, working_dir, input_file, dock_files_dir, output_dir, job_timeout_minutes=None):
        command_str = f"bash {self.SUBMISSION_SCRIPT_PATH}"
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

        qsub_args_str = f"-N {job_name} -o {output_dir.path} -e {output_dir.path} "

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
