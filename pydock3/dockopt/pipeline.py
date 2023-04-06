from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn, Iterable, Union
from datetime import datetime
import os
import functools
import re

import pandas as pd

from pydock3.files import Dir
from pydock3.dockopt.criterion import CRITERION_DICT
if TYPE_CHECKING:
    from pydock3.dockopt.results import ResultsManager


def add_timing_and_results_writing_to_run_method(pipeline_component: PipelineComponent) -> PipelineComponent:
    run = getattr(pipeline_component, "run")

    @functools.wraps(run)
    def new_run(
            self, 
            *args,
            components_to_run: str = "^.*$",  # default: match any string
            components_to_skip_if_results_exist: str = "^.*$",  # default: match any string
            components_to_force_rewrite_report: str = "^.*$",  # default: match any string,
            **kwargs
            ) -> pd.DataFrame:
        
        if re.match(components_to_run, rf"{self.component_id}") is None:
            return None

        if re.match(components_to_skip_if_results_exist, rf"{self.component_id}") is not None and self.results_manager is not None:
            if self.results_manager.results_exist(self):
                result = self.results_manager.load_results(self)
                if re.match(components_to_force_rewrite_report, rf"{self.component_id}") is not None:
                    self.results_manager.write_report(self)
                return result

        self.started_utc = datetime.utcnow()  # record utc datetime when `run` starts
        result = run(
            self, 
            *args, 
            components_to_run=components_to_run, 
            components_to_skip_if_results_exist=components_to_skip_if_results_exist, 
            components_to_force_rewrite_report=components_to_force_rewrite_report, 
            **kwargs
        )
        if self.results_manager is not None:  # write results if set
            self.results_manager.write_results(self, result)
            self.results_manager.write_report(self)
        self.finished_utc = datetime.utcnow()  # record utc datetime when `run` finishes

        return result

    setattr(pipeline_component, "run", new_run)

    return pipeline_component


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
        self.started_utc = None  # set by run()
        self.finished_utc = None  # set by run()

    def __init_subclass__(cls, **kwargs):
        return add_timing_and_results_writing_to_run_method(pipeline_component=cls)

    @property
    def component_dir(self):
        if self.component_id is None:  # component is pipeline
            return self.pipeline_dir
        else:
            return Dir(os.path.join(self.pipeline_dir.path, *self.component_id.split('.')))

    def run(
            self, 
            *args,
            components_to_run: str = "^.*$",  # default: match any string
            components_to_skip_if_results_exist: str = "^.*$",  # default: match any string
            components_to_force_rewrite_report: str = "^.*$",  # default: match any string,
            **kwargs
            ) -> NoReturn:
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
