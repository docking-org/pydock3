from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn, Iterable, Tuple
from datetime import datetime
import os

import pandas as pd

from pydock3.files import Dir
from pydock3.dockopt.criterion import criterion_dict
if TYPE_CHECKING:
    from pydock3.dockopt.results import ResultsManager


def add_timing_and_results_writing_to_run_method(_cls: object) -> object:
    run = getattr(_cls, "run")

    def new_run(self, *args, **kwargs):
        self.started_utc = datetime.utcnow()  # record utc datetime when `run` starts
        result = run(self, *args, **kwargs)
        if self.results_manager is not None:  # automatically write results if set
            self.results_manager.write_results(self, result)
        self.finished_utc = datetime.utcnow()  # record utc datetime when `run` finishes

        return result

    setattr(_cls, "run", new_run)

    return _cls


class PipelineObject(object):
    def __init__(
            self,
            job_dir_path: str,
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
    ):
        #
        self.job_dir = Dir(job_dir_path)
        self.top_n = top_n
        self.results_manager = results_manager

        #
        if criterion in criterion_dict:
            self.criterion = criterion_dict[criterion]()
        else:
            raise ValueError(f"`criterion` must be one of: {criterion_dict.keys()}. Witnessed: {criterion}")

        #
        self.started_utc = None  # set by add_timing_and_reporting_to_run_method() decorator; see .__init_subclass__()
        self.finished_utc = None  # set by add_timing_and_reporting_to_run_method() decorator; see .__init_subclass__()

    def __init_subclass__(cls, **kwargs):
        return add_timing_and_results_writing_to_run_method(_cls=cls)

    def run(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError

    def load_results_dataframe(self) -> pd.core.frame.DataFrame:
        return self.results_manager.load_results(self)

    def get_top_results_dataframe(self) -> pd.core.frame.DataFrame:
        return self.load_results_dataframe().nlargest(self.top_n, self.criterion.name)


class PipelineComponent(PipelineObject):
    def __init__(
            self,
            component_id: str,
            job_dir_path: str,
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
    ):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=results_manager,
        )

        #
        self.component_id = component_id

    def run(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError

    @property
    def component_dir(self):
        return Dir(os.path.join(self.job_dir.path, *self.component_id.split('.')))


class PipelineComponentSequenceIteration(PipelineComponent):
    def __init__(
        self,
        component_id: str,
        job_dir_path: str,
        criterion: str,
        top_n: int,
        results_manager: ResultsManager,
        components: Iterable[dict],
    ):
        super().__init__(
            component_id=component_id,
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=results_manager,
        )

        #
        self.components = components

    def run(self, *arg, **kwargs) -> NoReturn:
        raise NotImplementedError


class PipelineComponentSequence(PipelineComponent):
    def __init__(
        self,
        component_id: str,
        job_dir_path: str,
        criterion: str,
        top_n: int,
        results_manager: ResultsManager,
        components: Iterable[dict],
        num_repetitions: int,
        max_iterations_with_no_improvement: int,
        inter_iteration_criterion: str,
        inter_iteration_top_n: int,
    ):
        super().__init__(
            component_id=component_id,
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=results_manager,
        )

        #
        self.components = components
        self.num_repetitions = num_repetitions
        self.max_iterations_with_no_improvement = max_iterations_with_no_improvement
        self.inter_iteration_criterion = inter_iteration_criterion
        self.inter_iteration_top_n = inter_iteration_top_n

    def run(self, *arg, **kwargs) -> NoReturn:
        raise NotImplementedError


class Pipeline(PipelineObject):
    def __init__(
        self,
        job_dir_path: str,
        criterion: str,
        top_n: int,
        results_manager: ResultsManager,
        components: Iterable[dict],
    ):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=results_manager,
        )

        #
        self.components = components

    def run(self, *arg, **kwargs) -> NoReturn:
        raise NotImplementedError
