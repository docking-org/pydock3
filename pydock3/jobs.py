import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from enum import Enum

from pydock3.files import Dir, File, IndockFile
from pydock3.blastermaster.util import DockFiles
from pydock3.job_schedulers import JobScheduler

from pydock3.docking import __file__ as DOCKING_INIT_FILE_PATH

DOCK3_EXECUTABLE_PATH = os.path.join(
    os.path.dirname(DOCKING_INIT_FILE_PATH), "dock3", "dock64"
)
DOCK_RUN_SCRIPT_PATH = os.path.join(
    os.path.dirname(DOCKING_INIT_FILE_PATH), "rundock.bash"
)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
class JobSubmissionResult(Enum):
    SUCCESS = 1
    FAILED = 2
    SKIPPED_BECAUSE_ALREADY_COMPLETE = 3
    SKIPPED_BECAUSE_STILL_RUNNING = 4


@dataclass
class DockingJob(ABC):
    @abstractmethod
    def submit(self):
        raise NotImplementedError


@dataclass
class RetrodockJob(ABC):
    name: str
    job_dir: Dir
    input_sdi_file: File
    dock_files: DockFiles
    indock_file: IndockFile
    job_scheduler: JobScheduler
    temp_storage_path: str
    dock_executable_path: str = DOCK3_EXECUTABLE_PATH
    max_reattempts: int = 0
    num_attempts: int = 0

    N_TASKS = 2
    ACTIVES_TASK_ID = "1"
    DECOYS_TASK_ID = "2"
    OUTDOCK_FILE_NAME = "OUTDOCK.0"
    JOBLIST_FILE_NAME = "joblist"

    def __post_init__(self):
        #
        self.output_dir = Dir(
            path=os.path.join(self.job_dir.path, f"output"), create=True
        )

        #
        self.actives_output_dir = Dir(os.path.join(self.output_dir.path, self.ACTIVES_TASK_ID))
        self.decoys_output_dir = Dir(os.path.join(self.output_dir.path, self.DECOYS_TASK_ID))

        #
        self.actives_outdock_file = File(os.path.join(self.actives_output_dir.path, self.OUTDOCK_FILE_NAME))
        self.decoys_outdock_file = File(os.path.join(self.decoys_output_dir.path, self.OUTDOCK_FILE_NAME))

        #
        self._is_complete = False

    def submit(self, job_timeout_minutes=None, skip_if_complete=True):
        """
        if job submission is skipped, returns (JobSubmissionResult, None)
        if job submission is not skipped, returns (JobSubmissionResult, subprocess.CompletedProcess)
        """

        #
        if self.is_running:
            return JobSubmissionResult.SKIPPED_BECAUSE_STILL_RUNNING, None

        #
        if skip_if_complete:
            if self.is_complete:
                return JobSubmissionResult.SKIPPED_BECAUSE_ALREADY_COMPLETE, None

        # reset output dir
        self.output_dir.delete()
        self.output_dir.create()

        # make joblist
        with open(self.input_sdi_file.path, "r") as f:
            tgz_file_paths = [line.strip() for line in f.readlines()]
        try:
            assert len(tgz_file_paths) == 2
        except AssertionError:
            raise Exception(
                "Attempted to pass SDI file with more than two lines to RetrodockJob()"
            )
        with open(os.path.join(self.output_dir.path, self.JOBLIST_FILE_NAME), "w") as f:
            for i, tgz_file_path in enumerate(tgz_file_paths):
                task_id = i + 1
                f.write(f"{tgz_file_path} {task_id}\n")

        # set env vars dict
        dock_file_paths = [
            getattr(self.dock_files, dock_file_field.name).path
            for dock_file_field in fields(self.dock_files)
        ]
        env_vars_dict = {
            "EXPORT_DEST": self.output_dir.path,
            "INPUT_SOURCE": self.input_sdi_file.path,
            "DOCKEXEC": self.dock_executable_path,
            "TMPDIR": self.temp_storage_path,
            "DOCKFILE_PATHS_LIST": " ".join(dock_file_paths),
            "INDOCK_PATH": self.indock_file.path,
            "ONLY_EXPORT_MOL2_FOR_TASK_1": "true",
        }

        # submit job
        proc = self.job_scheduler.submit(
            job_name=self.name,
            script_path=DOCK_RUN_SCRIPT_PATH,
            env_vars_dict=env_vars_dict,
            output_dir_path=self.output_dir.path,
            n_tasks=self.N_TASKS,
            job_timeout_minutes=job_timeout_minutes,
        )
        self.num_attempts += 1

        if proc.stderr:
            return JobSubmissionResult.FAILED, proc
        else:
            return JobSubmissionResult.SUCCESS, proc

    @property
    def is_running(self):
        return self.job_scheduler.is_running_job(self.name)

    @property
    def is_complete(self):
        return all(
            [
                File.file_exists(
                    os.path.join(self.output_dir.path, task_id, self.OUTDOCK_FILE_NAME)
                )
                for task_id in [self.ACTIVES_TASK_ID, self.DECOYS_TASK_ID]
            ]
        )
