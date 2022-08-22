import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass

from pydock3.files import Dir, File
from pydock3.job_schedulers import JobScheduler

from pydock3.docking import __file__ as DOCKING_INIT_FILE_PATH
DOCK3_EXECUTABLE_PATH = os.path.join(os.path.dirname(DOCKING_INIT_FILE_PATH), 'dock3', 'dock64')


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


@dataclass
class DockingJob(ABC):

    @abstractmethod
    def run(self):
        raise NotImplementedError


@dataclass
class RetrospectiveDockingJob(ABC):
    name: str
    dock_files_dir: Dir
    input_file: File
    working_dir: Dir
    output_dir: Dir
    job_scheduler: JobScheduler
    dock_executable_path: str = DOCK3_EXECUTABLE_PATH
    max_reattempts: int = 0
    num_attempts: int = 0

    def __post_init__(self):
        self._is_complete = False

    def run(self, job_timeout_minutes=None, skip_if_complete=True):
        if skip_if_complete:
            if self.is_complete:
                return

        #
        self.output_dir.delete()
        self.output_dir.create()

        #
        proc = self.job_scheduler.run(
            job_name=self.name,
            dock_executable_path=self.dock_executable_path,
            working_dir=self.working_dir,
            input_file=self.input_file,
            dock_files_dir=self.dock_files_dir,
            output_dir=self.output_dir,
            job_timeout_minutes=job_timeout_minutes,
        )
        self.num_attempts += 1

        return proc

    @property
    def is_running(self):
        return self.job_scheduler.is_running_job(self.name)

    @property
    def is_complete(self):
        return all([File.file_exists(os.path.join(self.output_dir.path, sub_job_num, "OUTDOCK.0")) for sub_job_num in ['1', '2']])
