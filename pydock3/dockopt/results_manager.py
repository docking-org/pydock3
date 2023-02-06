from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os

import pandas as pd

from pydock3.dockopt.reporter import PDFReporter
if TYPE_CHECKING:
    from pydock3.dockopt.pipeline import PipelineComponent


#
RESULTS_CSV_FILE_NAME = "results.csv"


def add_sorting_of_results_dataframe_by_criterion_to_write_results_method(_cls: object) -> object:
    write_results = getattr(_cls, "write_results")

    def new_write_results(self, pipeline_component):
        return write_results(self, pipeline_component).sort_values(by=pipeline_component.criterion.name, ascending=False, ignore_index=True)

    setattr(_cls, "write_results", new_write_results)

    return _cls


class ResultsManager(object):
    def __init__(self, results_file_name: str):
        self.results_file_name = results_file_name

    def __init_subclass__(cls, **kwargs):
        return add_sorting_of_results_dataframe_by_criterion_to_write_results_method(_cls=cls)

    def write_results(
        self,
        pipeline_component: PipelineComponent,
        results_dataframe: pd.core.frame.DataFrame,
    ) -> NoReturn:
        raise NotImplementedError

    def load_results(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError

    def write_report(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError


class DockoptResultsManager(ResultsManager):
    def __init__(self, results_file_name: str):
        super().__init__(results_file_name)

    def write_results(
        self,
        pipeline_component: PipelineComponent,
        results_dataframe: pd.core.frame.DataFrame,
    ) -> None:
        results_dataframe.to_csv(os.path.join(pipeline_component.dir.path, self.results_file_name))

    def load_results(self, pipeline_component: PipelineComponent) -> pd.core.frame.DataFrame:
        df = pd.read_csv(os.path.join(pipeline_component.dir.path, self.results_file_name))
        df = df.loc[
            :, ~df.columns.str.contains("^Unnamed")
        ]  # remove useless index column

        return df

    def write_report(self, pipeline_component: PipelineComponent) -> None:
        PDFReporter().write_report(pipeline_component)
