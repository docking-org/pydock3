import logging
from typing import Union, List, Iterable
import os
from abc import ABC, abstractmethod
from itertools import groupby
from operator import itemgetter
import re
from subprocess import CompletedProcess

import xmltodict

from pydock3.util import system_call, get_nested_dict_item
from pydock3.files import File

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
    def job_is_on_queue(self, job_name: str) -> bool:
        raise NotImplementedError

    @abstractmethod
    def task_is_on_queue(self, task_id: Union[str, int], job_name: str) -> bool:
        raise NotImplementedError


class SlurmJobScheduler(JobScheduler):
    REQUIRED_ENV_VAR_NAMES = [
        "SBATCH_EXEC",
        "SQUEUE_EXEC",
    ]

    MAX_ARRAY_JOBS_SUBMITTED_PER_SBATCH = 100

    def __init__(self) -> None:
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
            job_name: str,
            script_path: str,
            env_vars_dict: dict,
            out_log_dir_path: str,
            err_log_dir_path: str,
            task_ids: Iterable[Union[str, int]],
            job_timeout_minutes: Union[int, None] = None,
    ) -> List[CompletedProcess]:
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
            if self.SLURM_SETTINGS:
                if File.file_exists(self.SLURM_SETTINGS):  # TODO: move validation to __init__
                    command_str = f"source {self.SLURM_SETTINGS}; {command_str}"
            if job_timeout_minutes is None:
                command_str += " --time=0"
            else:
                command_str += f" --time={job_timeout_minutes}"
            proc = system_call(
                command_str, env_vars_dict=env_vars_dict
            )  # need to pass env_vars_dict here so that '--export=ALL' in command can pass along all the env vars
            procs.append(proc)

        return procs

    def job_is_on_queue(self, job_name: str) -> bool:
        command_str = f"{self.SQUEUE_EXEC} --format='%i %j %t' | grep '{job_name}'"
        proc = system_call(command_str)

        #
        if proc.stdout:
            return True
        else:
            return False

    def task_is_on_queue(self, task_id: Union[str, int], job_name: str) -> bool:
        command_str = f"{self.SQUEUE_EXEC} -r --format='%i %j %t' | grep '{job_name}'"
        proc = system_call(command_str)

        #
        if not proc.stdout:
            return False

        #
        for line in proc.stdout.split('\n'):
            line_stripped = line.strip()
            if line_stripped:
                job_id, job_name, state = line_stripped.split()
                if job_id.endswith(f"_{task_id}"):
                    return True

        return False


class SGEJobScheduler(JobScheduler):
    REQUIRED_ENV_VAR_NAMES = [
        "QSUB_EXEC",
        "QSTAT_EXEC",
    ]

    def __init__(self) -> None:
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
            job_name: str,
            script_path: str,
            env_vars_dict: dict,
            out_log_dir_path: str,
            err_log_dir_path: str,
            task_ids: Iterable[Union[str, int]],
            job_timeout_minutes: Union[int, None] = None,
    ) -> List[CompletedProcess]:
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
            command_str = f"{self.QSUB_EXEC} -V -N {job_name} -o {out_log_dir_path} -e {err_log_dir_path} -cwd -S /bin/bash -q !gpu.q -t {array_str} {script_path}"
            if self.SGE_SETTINGS:
                if File.file_exists(self.SGE_SETTINGS):  # TODO: move validation to __init__
                    command_str = f"source {self.SGE_SETTINGS}; {command_str}"
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

    def job_is_on_queue(self, job_name: str) -> bool:
        command_str = f"{self.QSTAT_EXEC} -r | grep '{job_name}'"
        proc = system_call(command_str)
        if proc.stdout:
            return True
        else:
            return False

    def _get_qstat_xml_as_dict(self) -> dict:
        command_str = f"{self.QSTAT_EXEC} -xml"
        proc = system_call(command_str)
        if proc.stdout is not None:
            return xmltodict.parse(proc.stdout)
        else:
            raise Exception(f"Command '{command_str}' returned stdout of None. stderr: {proc.stderr}")

    def task_is_on_queue(self, task_id: Union[str, int], job_name: str) -> bool:
        task_num = int(task_id)

        #
        q_dict = self._get_qstat_xml_as_dict()

        #
        try:
            obj = get_nested_dict_item(q_dict, ('job_info', 'queue_info', 'job_list'))
        except (KeyError, TypeError):
            return False

        #
        if isinstance(obj, dict):
            job_dicts = [obj]
        elif isinstance(obj, list):
            job_dicts = obj
        else:
            raise Exception(f"Unexpected type for 'job_list': {type(obj)}")

        #
        for job_dict in job_dicts:
            #
            if not job_dict.get('JB_name') == job_name:
                continue

            #
            tasks_str = job_dict.get('tasks')
            if tasks_str is None:
                continue

            #
            pattern = r'^(\d+)-(\d+)(:\d+)?$'
            match = re.match(pattern, tasks_str)
            if match is not None:
                if match.group(1) is not None and match.group(2) is not None:
                    start = int(match.group(1))
                    end = int(match.group(2))
                    if (task_num >= start) and (task_num <= end):
                        #
                        return True

            #
            pattern = r'^(\d+)$'
            match = re.match(pattern, tasks_str)
            if match is not None:
                if match.group(1) is not None:
                    num = int(match.group(1))
                    if task_num == num:
                        return True

        #
        return False
