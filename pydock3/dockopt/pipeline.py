from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn, Iterable, Tuple
from datetime import datetime

import pandas as pd

from pydock3.files import Dir
from pydock3.dockopt.criterion import criterion_dict
if TYPE_CHECKING:
    from pydock3.dockopt.results_manager import ResultsManager


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


class PipelineComponent(object):
    def __init__(
            self,
            component_id: str,
            dir_path: str,
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
    ):
        #
        self.component_id = component_id
        self.dir = Dir(dir_path)
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

        #
        # TODO: create directory

    def __init_subclass__(cls, **kwargs):
        return add_timing_and_results_writing_to_run_method(_cls=cls)

    def run(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError

    def load_results_dataframe(self) -> pd.core.frame.DataFrame:
        return self.results_manager.load_results(self)

    def get_top_results_dataframe(self) -> pd.core.frame.DataFrame:
        return self.load_results_dataframe().nlargest(self.top_n, self.criterion.name)


class PipelineComponentSequenceIteration(PipelineComponent):
    def __init__(
        self,
        component_id: str,
        dir_path: str,
        criterion: str,
        top_n: int,
        results_manager: ResultsManager,
        components: Iterable[dict],
    ):
        super().__init__(
            component_id=component_id,
            dir_path=dir_path,
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
        dir_path: str,
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
            dir_path=dir_path,
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
