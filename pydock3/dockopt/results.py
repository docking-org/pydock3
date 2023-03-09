from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os
import shutil
import glob
import logging
from dataclasses import fields

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joypy import joyplot

from pydock3.files import Dir
from pydock3.dockopt.util import BEST_RETRODOCK_JOBS_DIR_NAME
from pydock3.dockopt.docking_configuration import DockingConfiguration
from pydock3.dockopt.roc import ROC
from pydock3.dockopt.reporter import PDFReporter
from pydock3.retrodock.retrodock import ROC_PLOT_FILE_NAME, ENERGY_PLOT_FILE_NAME, CHARGE_PLOT_FILE_NAME, str_to_float, get_results_dataframe_from_actives_job_and_decoys_job_outdock_files
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
        self.write_report(pipeline_component)

    def load_results(self, pipeline_component: PipelineComponent) -> pd.DataFrame:
        df = pd.read_csv(os.path.join(pipeline_component.component_dir.path, self.results_file_name))
        df = df.loc[
            :, ~df.columns.str.contains("^Unnamed")
        ]  # remove useless index column

        return df

    def write_report(self, pipeline_component: PipelineComponent) -> None:
        PDFReporter().write_report(pipeline_component)

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

            # copy docking configuration files to best jobs dir
            dst_best_job_dir_path = os.path.join(pipeline_component.best_retrodock_jobs_dir.path, f"rank={i+1}_step={dc.component_id}_conf={dc.configuration_num}")
            best_job_dockfiles_dir = Dir(
                os.path.join(dst_best_job_dir_path, "dockfiles"),
                create=True
            )
            dock_files = dc.get_dock_files(pipeline_component.pipeline_dir.path)
            for field in fields(dock_files):
                best_job_dockfiles_dir.copy_in_file(getattr(dock_files, field.name).path)
            best_job_dockfiles_dir.copy_in_file(dc.get_indock_file(pipeline_component.pipeline_dir.path).path)

            #
            src_retrodock_job_actives_dir_path = os.path.join(pipeline_component.retrodock_jobs_dir.path, "actives", str(dc.configuration_num))
            src_retrodock_job_decoys_dir_path = os.path.join(pipeline_component.retrodock_jobs_dir.path, "decoys", str(dc.configuration_num))
            shutil.copytree(
                src_retrodock_job_actives_dir_path,
                os.path.join(dst_best_job_dir_path, "actives"),
            )
            shutil.copytree(
                src_retrodock_job_decoys_dir_path,
                os.path.join(dst_best_job_dir_path, "decoys"),
            )

            #
            df_best_job = (
                get_results_dataframe_from_actives_job_and_decoys_job_outdock_files(
                    actives_outdock_file_path=os.path.join(
                        dst_best_job_dir_path,
                        "actives",
                        "OUTDOCK.0",
                    ),
                    decoys_outdock_file_path=os.path.join(
                        dst_best_job_dir_path,
                        "decoys",
                        "OUTDOCK.0",
                    ),
                )
            )

            # sort dataframe by total energy score
            df_best_job = df_best_job.sort_values(
                by=["total_energy", "is_active"],
                na_position="last", ignore_index=True
            )  # sorting secondarily by 'is_active' (0 or 1) ensures that decoys are ranked before actives in case they have the same exact score (pessimistic approach)
            df_best_job = df_best_job.drop_duplicates(
                subset=["db2_file_path"], keep="first",
                ignore_index=True
            )  # keep only the best score per molecule

            # get ROC and calculate enrichment score of this job's docking set-up
            # TODO: get this from Retrodock instead
            booleans = df_best_job["is_active"].astype(bool)
            roc = ROC(booleans)

            # write ROC plot image
            roc_plot_image_path = os.path.join(dst_best_job_dir_path, ROC_PLOT_FILE_NAME)
            roc.plot(save_path=roc_plot_image_path)

            # ridge plot for energy terms
            # TODO: get this from Retrodock instead
            pivot_rows = []
            for i, best_job_row in df_best_job.iterrows():
                for col in [
                    "total_energy",
                    "electrostatic_energy",
                    "vdw_energy",
                    "polar_desolvation_energy",
                    "apolar_desolvation_energy",
                ]:
                    pivot_row = {"energy_term": col}
                    if best_job_row["is_active"] == 1:
                        pivot_row["active"] = str_to_float(best_job_row[col])
                        pivot_row["decoy"] = np.nan
                    else:
                        pivot_row["active"] = np.nan
                        pivot_row["decoy"] = str_to_float(best_job_row[col])
                    pivot_rows.append(pivot_row)
            df_best_job_pivot = pd.DataFrame(pivot_rows)
            fig, ax = joyplot(
                data=df_best_job_pivot,
                by="energy_term",
                column=["active", "decoy"],
                color=["#686de0", "#eb4d4b"],
                legend=True,
                alpha=0.85,
                figsize=(12, 8),
                ylim="own",
            )
            plt.title("ridgeline plot: energy terms (actives vs. decoys)")
            plt.tight_layout()
            plt.savefig(os.path.join(dst_best_job_dir_path, ENERGY_PLOT_FILE_NAME))
            plt.close(fig)

            # split violin plot of charges
            # TODO: get this from Retrodock instead
            fig = plt.figure()
            sns.violinplot(
                data=df_best_job,
                x="charge",
                y="total_energy",
                split=True,
                hue="activity_class",
            )
            plt.title("split violin plot: charge (actives vs. decoys)")
            plt.tight_layout()
            plt.savefig(os.path.join(dst_best_job_dir_path, CHARGE_PLOT_FILE_NAME))
            plt.close(fig)


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
            src_best_job_dir_path, = tuple(glob.glob(os.path.join(pipeline_component.pipeline_dir.path, *dc.component_id.split('.'), BEST_RETRODOCK_JOBS_DIR_NAME, f"rank=*_step={dc.component_id}_conf={dc.configuration_num}")))
            dst_best_job_dir_path = os.path.join(
                pipeline_component.best_retrodock_jobs_dir.path,
                f"rank={i + 1}_step={dc.component_id}_conf={dc.configuration_num}",
            )
            shutil.copytree(
                src_best_job_dir_path,
                dst_best_job_dir_path,
            )


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
            src_best_job_dir_path, = tuple(glob.glob(os.path.join(pipeline_component.pipeline_dir.path, *dc.component_id.split('.'), BEST_RETRODOCK_JOBS_DIR_NAME, f"rank=*_step={dc.component_id}_conf={dc.configuration_num}")))
            dst_best_job_dir_path = os.path.join(
                pipeline_component.best_retrodock_jobs_dir.path,
                f"rank={i + 1}_step={dc.component_id}_conf={dc.configuration_num}",
            )
            shutil.copytree(
                src_best_job_dir_path,
                dst_best_job_dir_path,
            )
