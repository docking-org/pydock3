import logging
import os
from abc import ABC, abstractmethod
from itertools import groupby
from operator import itemgetter

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
            out_log_dir_path,
            err_log_dir_path,
            task_ids,
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

    MAX_ARRAY_JOBS_SUBMITTED_PER_SBATCH = 100

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
            out_log_dir_path,
            err_log_dir_path,
            task_ids,
            job_timeout_minutes=None,
    ):
        #
        task_nums = sorted([int(task_id) for task_id in task_ids])
        contiguous_task_nums_sets = [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(task_nums), lambda x: x[0] - x[1])]

        #
        procs = []
        for contiguous_task_nums_set in contiguous_task_nums_sets:
            if len(contiguous_task_nums_set) == 1:
                array_str = f"{contiguous_task_nums_set[0]}"
            else:
                array_str = f"{contiguous_task_nums_set[0]}-{contiguous_task_nums_set[-1]}"
            command_str = f"{self.SBATCH_EXEC} --export=ALL -J {job_name} -o {out_log_dir_path}/{job_name}_%A_%a.out -e {err_log_dir_path}/{job_name}_%A_%a.err --signal=B:USR1@120 --array={array_str} {script_path}"
            if job_timeout_minutes is None:
                command_str += " --time=0"
            else:
                command_str += f" --time={job_timeout_minutes}"
            proc = system_call(
                command_str, env_vars_dict=env_vars_dict
            )  # need to pass env_vars_dict here so that '--export=ALL' in command can pass along all the env vars
            procs.append(proc)

        return procs

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
            out_log_dir_path,
            err_log_dir_path,
            task_ids,
            job_timeout_minutes=None,
    ):
        #
        if not job_name[0].isalpha():
            raise Exception(f"{self.name} job names must start with a letter.")

        #
        task_nums = sorted([int(task_id) for task_id in task_ids])
        contiguous_task_nums_sets = [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(task_nums), lambda x: x[0] - x[1])]

        #
        procs = []
        for contiguous_task_nums_set in contiguous_task_nums_sets:
            if len(contiguous_task_nums_set) == 1:
                array_str = f"{contiguous_task_nums_set[0]}"
            else:
                array_str = f"{contiguous_task_nums_set[0]}-{contiguous_task_nums_set[-1]}"
            command_str = f"source {self.SGE_SETTINGS}; {self.QSUB_EXEC} -V -N {job_name} -o {out_log_dir_path} -e {err_log_dir_path} -cwd -S /bin/bash -q !gpu.q -t {array_str} {script_path}"
            if job_timeout_minutes is not None:
                job_timeout_seconds = 60 * job_timeout_minutes
                command_str += (
                    f" -l s_rt={job_timeout_seconds} -l h_rt={job_timeout_seconds} "
                )
            proc = system_call(
                command_str, env_vars_dict=env_vars_dict
            )  # need to pass env_vars_dict here so that '-V' in command can pass along all the env vars
            procs.append(proc)

        return procs

    def is_running_job(self, job_name):
        command_str = f"{self.QSTAT_EXEC} -r | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False
