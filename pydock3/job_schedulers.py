import logging
import os
from abc import ABC, abstractmethod

from pydock3.util import system_call

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class JobScheduler(ABC):
    REQUIRED_ENV_VAR_NAMES = []

    def __init__(self, name):
        self.name = name

    @abstractmethod
    def submit(
            self,
            job_name,
            script_path,
            env_vars_dict,
            output_dir_path,
            n_tasks,
            job_timeout_minutes=None,
    ):
        """returns: subprocess.CompletedProcess"""

        raise NotImplementedError

    @abstractmethod
    def is_running_job(self, job_name):
        raise NotImplementedError



class SlurmJobScheduler(JobScheduler):
    REQUIRED_ENV_VAR_NAMES = [
        "SBATCH_EXEC",
        "SQUEUE_EXEC",
    ]

    def __init__(self):
        super().__init__(name="Slurm")

        # set required env vars
        self.SBATCH_EXEC = os.environ["SBATCH_EXEC"]
        self.SQUEUE_EXEC = os.environ["SQUEUE_EXEC"]

        # set optional env vars
        self.SLURM_SETTINGS = os.environ.get("SLURM_SETTINGS")

        #
        if self.SLURM_SETTINGS:
            system_call(f"source {self.SLURM_SETTINGS}")

    def submit(
            self,
            job_name,
            script_path,
            env_vars_dict,
            output_dir_path,
            n_tasks,
            job_timeout_minutes=None,
    ):
        command_str = f"{self.SBATCH_EXEC} --export=ALL -J {job_name} -o {output_dir_path}/{job_name}_%A_%a.out -e {output_dir_path}/{job_name}_%A_%a.err --signal=B:USR1@120 --array=1-{n_tasks} {script_path}"
        if job_timeout_minutes is None:
            command_str += " --time=0"
        else:
            command_str += f" --time={job_timeout_minutes}"

        return system_call(
            command_str, env_vars_dict=env_vars_dict
        )  # need to pass env_vars_dict here so that '--export=ALL' in command can pass along all the env vars

    def is_running_job(self, job_name):
        command_str = f"{self.SQUEUE_EXEC} --format='%.18i %.{len(job_name)}j' | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False


class SGEJobScheduler(JobScheduler):
    REQUIRED_ENV_VAR_NAMES = [
        "QSUB_EXEC",
        "QSTAT_EXEC",
    ]

    def __init__(self):
        super().__init__(name="SGE")

        # set required env vars
        self.QSUB_EXEC = os.environ["QSUB_EXEC"]
        self.QSTAT_EXEC = os.environ["QSTAT_EXEC"]

        # set optional env vars
        self.SGE_SETTINGS = os.environ.get("SGE_SETTINGS")

        #
        if self.SGE_SETTINGS:
            system_call(f"source {self.SGE_SETTINGS}")

    def submit(
            self,
            job_name,
            script_path,
            env_vars_dict,
            output_dir_path,
            n_tasks,
            job_timeout_minutes=None,
    ):
        #
        if not job_name[0].isalpha():
            raise Exception(f"{self.name} job names must start with a letter.")

        #
        command_str = f"source {self.SGE_SETTINGS}; {self.QSUB_EXEC} -V -N {job_name} -o {output_dir_path} -e {output_dir_path} -cwd -S /bin/bash -q !gpu.q -t 1-{n_tasks} {script_path}"
        if job_timeout_minutes is not None:
            job_timeout_seconds = 60 * job_timeout_minutes
            command_str += (
                f" -l s_rt={job_timeout_seconds} -l h_rt={job_timeout_seconds} "
            )

        return system_call(
            command_str, env_vars_dict=env_vars_dict
        )  # need to pass env_vars_dict here so that '-V' in command can pass along all the env vars

    def is_running_job(self, job_name):
        command_str = f"{self.QSTAT_EXEC} -r | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False
