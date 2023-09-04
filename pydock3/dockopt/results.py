from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os
import glob
import logging
from dataclasses import fields

import pandas as pd

from pydock3.files import Dir, create_relative_symlink
from pydock3.dockopt.util import BEST_RETRODOCK_JOBS_DIR_NAME
from pydock3.dockopt.docking_configuration import DockingConfiguration
from pydock3.jobs import OUTDOCK_FILE_NAME
from pydock3.dockopt.reporter import HTMLReporter
from pydock3.retrodock.retrodock import ROC_PLOT_FILE_NAME, ENERGY_TERMS_PLOT_FILE_NAME, CHARGE_PLOT_FILE_NAME, str_to_float, get_results_dataframe_from_positives_job_and_negatives_job_outdock_files, process_retrodock_job_results
if TYPE_CHECKING:
    from pydock3.dockopt.pipeline import PipelineComponent


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def add_sorting_of_results_dataframe_by_criterion_to_write_results_method(_cls: object) -> object:
    write_results = getattr(_cls, "write_results")

    def new_write_results(self, pipeline_component, results_dataframe):
        return write_results(self, pipeline_component, results_dataframe.sort_values(by=pipeline_component.criterion.name, ascending=False, ignore_index=True))

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
        results_dataframe: pd.DataFrame,
    ) -> NoReturn:
        raise NotImplementedError
    
    def results_exist(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError

    def load_results(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError

    def write_report(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError


class DockoptPipelineComponentResultsManager(ResultsManager):
    def __init__(self, results_file_name: str):
        super().__init__(results_file_name)

    def write_results(
        self,
        pipeline_component: PipelineComponent,
        results_dataframe: pd.DataFrame,
    ) -> None:
        results_dataframe.to_csv(os.path.join(pipeline_component.component_dir.path, self.results_file_name))
        self.save_best_retrodock_jobs(pipeline_component)

    def results_exist(self, pipeline_component: PipelineComponent) -> bool:
        return os.path.exists(os.path.join(pipeline_component.component_dir.path, self.results_file_name))

    def load_results(self, pipeline_component: PipelineComponent) -> pd.DataFrame:
        df = pd.read_csv(os.path.join(pipeline_component.component_dir.path, self.results_file_name))
        df = df.loc[
            :, ~df.columns.str.contains("^Unnamed")
        ]  # remove useless index column

        return df

    def write_report(self, pipeline_component: PipelineComponent) -> None:
        HTMLReporter().write_report(pipeline_component)

    def save_best_retrodock_jobs(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError


class DockoptStepResultsManager(DockoptPipelineComponentResultsManager):

    def __init__(self, results_file_name: str):
        super().__init__(results_file_name)

    def save_best_retrodock_jobs(self, pipeline_component: PipelineComponent):
        # reset best jobs dir
        pipeline_component.best_retrodock_jobs_dir.reset()

        #
        logger.debug(
            f"Copying top {pipeline_component.top_n} retrodock jobs to {pipeline_component.best_retrodock_jobs_dir.path}"
        )
        for i, row in pipeline_component.get_top_results_dataframe().iterrows():
            #
            dc = DockingConfiguration.from_dict(row.to_dict())

            #
            dst_best_job_dir = Dir(os.path.join(pipeline_component.best_retrodock_jobs_dir.path, f"rank={i+1}-step={dc.component_id}-conf={dc.configuration_num}"), create=True, reset=True)
            best_job_dockfiles_dir = Dir(
                os.path.join(dst_best_job_dir.path, "dockfiles"),
                create=True,
                reset=True,
            )  # create best job dockfiles dir

            # create symbolic links instead of copying in order to save time & space
            dock_files = dc.get_dock_files(pipeline_component.pipeline_dir.path)
            for field in fields(dock_files):
                dock_file = getattr(dock_files, field.name)
                create_relative_symlink(dock_file.path, os.path.join(best_job_dockfiles_dir.path, dock_file.name), target_is_directory=False)

            indock_file = dc.get_indock_file(pipeline_component.pipeline_dir.path)
            create_relative_symlink(indock_file.path, os.path.join(best_job_dockfiles_dir.path, indock_file.name), target_is_directory=False)

            src_retrodock_job_positives_dir_path = os.path.join(pipeline_component.retrodock_jobs_dir.path, "positives", str(dc.configuration_num))
            src_retrodock_job_negatives_dir_path = os.path.join(pipeline_component.retrodock_jobs_dir.path, "negatives", str(dc.configuration_num))

            dst_retrodock_job_positives_dir_path = os.path.join(dst_best_job_dir.path, "positives")
            dst_retrodock_job_negatives_dir_path = os.path.join(dst_best_job_dir.path, "negatives")

            create_relative_symlink(
                src_retrodock_job_positives_dir_path,
                dst_retrodock_job_positives_dir_path,
                target_is_directory=True,
            )
            create_relative_symlink(
                src_retrodock_job_negatives_dir_path,
                dst_retrodock_job_negatives_dir_path,
                target_is_directory=True,
            )

            # save plots and other info
            process_retrodock_job_results(
                positives_outdock_file_path=os.path.join(dst_retrodock_job_positives_dir_path, OUTDOCK_FILE_NAME),
                negatives_outdock_file_path=os.path.join(dst_retrodock_job_negatives_dir_path, OUTDOCK_FILE_NAME),
                save_dir_path=dst_best_job_dir.path,
            )


class DockoptStepSequenceIterationResultsManager(DockoptPipelineComponentResultsManager):

    def __init__(self, results_file_name: str):
        super().__init__(results_file_name)

    def save_best_retrodock_jobs(self, pipeline_component: PipelineComponent):
        # reset best jobs dir
        pipeline_component.best_retrodock_jobs_dir.reset()

        #
        logger.debug(
            f"Copying top {pipeline_component.top_n} retrodock jobs to {pipeline_component.best_retrodock_jobs_dir.path}"
        )
        for i, row in pipeline_component.get_top_results_dataframe().iterrows():
            #
            dc = DockingConfiguration.from_dict(row.to_dict())

            #
            src_best_job_dir_path, = tuple(glob.glob(os.path.join(pipeline_component.pipeline_dir.path, *dc.component_id.split('.'), BEST_RETRODOCK_JOBS_DIR_NAME, f"rank=*-step={dc.component_id}-conf={dc.configuration_num}")))
            dst_best_job_dir_path = os.path.join(
                pipeline_component.best_retrodock_jobs_dir.path,
                f"rank={i + 1}-step={dc.component_id}-conf={dc.configuration_num}",
            )

            # create a symbolic link instead of copying in order to save time & space
            create_relative_symlink(src_best_job_dir_path, dst_best_job_dir_path, target_is_directory=True)


class DockoptStepSequenceResultsManager(DockoptPipelineComponentResultsManager):

    def __init__(self, results_file_name: str):
        super().__init__(results_file_name)

    def save_best_retrodock_jobs(self, pipeline_component: PipelineComponent):
        # reset best jobs dir
        pipeline_component.best_retrodock_jobs_dir.reset()

        #
        logger.debug(
            f"Copying top {pipeline_component.top_n} retrodock jobs to {pipeline_component.best_retrodock_jobs_dir.path}"
        )
        for i, row in pipeline_component.get_top_results_dataframe().iterrows():
            #
            dc = DockingConfiguration.from_dict(row.to_dict())

            #
            src_best_job_dir_path, = tuple(glob.glob(os.path.join(pipeline_component.pipeline_dir.path, *dc.component_id.split('.'), BEST_RETRODOCK_JOBS_DIR_NAME, f"rank=*-step={dc.component_id}-conf={dc.configuration_num}")))
            dst_best_job_dir_path = os.path.join(
                pipeline_component.best_retrodock_jobs_dir.path,
                f"rank={i + 1}-step={dc.component_id}-conf={dc.configuration_num}",
            )

            # create a symbolic link instead of copying in order to save time & space
            create_relative_symlink(src_best_job_dir_path, dst_best_job_dir_path, target_is_directory=True)
