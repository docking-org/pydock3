from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn, Iterable, Union
from datetime import datetime
import os

import pandas as pd

from pydock3.files import Dir
from pydock3.dockopt.criterion import CRITERION_DICT
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


class PipelineComponent(object):
    def __init__(
            self,
            pipeline_dir_path: str,
            component_id: Union[str, None],
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
    ):
        #
        self.pipeline_dir = Dir(pipeline_dir_path)
        self.component_id = component_id
        self.top_n = top_n
        self.results_manager = results_manager

        #
        if criterion in CRITERION_DICT:
            self.criterion = CRITERION_DICT[criterion]()
        else:
            raise ValueError(f"`criterion` must be one of: {CRITERION_DICT.keys()}. Witnessed: {criterion}")

        #
        self.started_utc = None  # set by add_timing_and_reporting_to_run_method() decorator; see .__init_subclass__()
        self.finished_utc = None  # set by add_timing_and_reporting_to_run_method() decorator; see .__init_subclass__()

    def __init_subclass__(cls, **kwargs):
        return add_timing_and_results_writing_to_run_method(_cls=cls)

    @property
    def component_dir(self):
        if self.component_id is None:  # component is pipeline
            return self.pipeline_dir
        else:
            return Dir(os.path.join(self.pipeline_dir.path, *self.component_id.split('.')))

    def run(self, *args, **kwargs) -> NoReturn:
        raise NotImplementedError

    def load_results_dataframe(self) -> pd.DataFrame:
        return self.results_manager.load_results(self)

    def get_top_results_dataframe(self) -> pd.DataFrame:
        return self.load_results_dataframe().nlargest(self.top_n, self.criterion.name)


class PipelineComponentSequenceIteration(PipelineComponent):
    def __init__(
            self,
            pipeline_dir_path: str,
            component_id: str,
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
            components: Iterable[dict],
    ):
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
            component_id=component_id,
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
            pipeline_dir_path: str,
            component_id: str,
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
            pipeline_dir_path=pipeline_dir_path,
            component_id=component_id,
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


class Pipeline(PipelineComponent):
    def __init__(
            self,
            pipeline_dir_path: str,
            criterion: str,
            top_n: int,
            results_manager: ResultsManager,
            components: Iterable[dict],
    ):
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
            component_id=None,
            criterion=criterion,
            top_n=top_n,
            results_manager=results_manager,
        )

        #
        self.components = components

    def run(self, *arg, **kwargs) -> NoReturn:
        raise NotImplementedError
