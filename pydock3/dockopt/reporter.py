from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os
import itertools

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from pydock3.files import File
from pydock3.retrodock.retrodock import ROC_PLOT_FILE_NAME, ENERGY_PLOT_FILE_NAME, CHARGE_PLOT_FILE_NAME
from pydock3.dockopt.criterion import criterion_dict
if TYPE_CHECKING:
    from pydock3.dockopt.pipeline import PipelineComponent


#
RETRODOCK_JOB_DIR_PATH_COLUMN_NAME = "retrodock_job_dir"

#
METRICS = list(criterion_dict.keys())
ALL_POSSIBLE_NON_PARAMETER_COLUMNS = [RETRODOCK_JOB_DIR_PATH_COLUMN_NAME] + METRICS


class Reporter(object):
    def __init__(self, report_file_name: str):
        self.report_file_name = report_file_name

    def write_report(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError


class PDFReporter(Reporter):
    def __init__(self, report_file_name: str = "report.pdf"):
        super().__init__(report_file_name)

        #
        plt.rcParams.update({"font.size": 14})

    def write_report(self, pipeline_component: PipelineComponent) -> None:
        def get_ordinal(n):
            return "%d%s" % (
                n,
                "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10:: 4],
            )

        df = pipeline_component.load_results_dataframe()

        with PdfPages(os.path.join(pipeline_component.dir.path, self.report_file_name)) as f:
            #
            fig = plt.figure(figsize=(11.0, 8.5))

            # write box plot for data grouped by runnable sequence identifier
            if "pipeline_component_id" in df.columns:
                # TODO
                pass

            #
            for i, best_job_dir_path in enumerate(
                    pipeline_component.get_top_results_dataframe()):  # TODO: replace
                #
                roc_file_path = os.path.join(
                    best_job_dir_path, ROC_PLOT_FILE_NAME
                )
                if File.file_exists(
                        roc_file_path):  # TODO: replace this logic with call to `Retrodock.analyze` method or something in order to ensure that these Retrodock-pertaining files are created
                    image = mpimg.imread(roc_file_path)
                    plt.axis("off")
                    plt.suptitle(
                        f"linear-log ROC plot of {get_ordinal(i + 1)} best job\n{best_job_dir_path}"
                    )
                    plt.imshow(image)
                    f.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                #
                energy_plot_file_path = os.path.join(
                    best_job_dir_path, ENERGY_PLOT_FILE_NAME
                )
                if File.file_exists(energy_plot_file_path):
                    fig = plt.figure(figsize=(11.0, 8.5))
                    image = mpimg.imread(energy_plot_file_path)
                    plt.axis("off")
                    plt.suptitle(
                        "energy terms ridgeline plot of best job")
                    plt.imshow(image)
                    f.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                #
                charge_plot_file_path = os.path.join(
                    best_job_dir_path, CHARGE_PLOT_FILE_NAME
                )
                if File.file_exists(charge_plot_file_path):
                    fig = plt.figure(figsize=(11.0, 8.5))
                    image = mpimg.imread(charge_plot_file_path)
                    plt.axis("off")
                    plt.suptitle("charge violin plot of best job")
                    plt.imshow(image)
                    f.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

            #
            multivalued_config_param_columns = [
                column
                for column in df.columns
                if column not in ALL_POSSIBLE_NON_PARAMETER_COLUMNS
                   and df[column].nunique() > 1
            ]
            for column in multivalued_config_param_columns:
                fig = plt.figure(figsize=(8, 6))

                if len(df) == len(
                        df.drop_duplicates(subset=[column])):
                    df.plot.bar(x=column, y=pipeline_component.criterion.name,
                                rot=0)
                else:
                    sns.boxplot(
                        data=df,
                        x=column,
                        y=pipeline_component.criterion.name,
                        showfliers=False,
                        boxprops={"facecolor": "None"},
                    )
                    sns.stripplot(data=df, x=column,
                                  y=pipeline_component.criterion.name, zorder=0.5)
                    fig.autofmt_xdate(rotation=25)
                    plt.yticks(rotation=0)

                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)

            #
            df = df.sort_values(
                by=pipeline_component.criterion.name, ascending=False,
                ignore_index=True
            )
            for column_1, column_2 in itertools.combinations(
                    multivalued_config_param_columns, 2
            ):
                #
                fig, ax = plt.subplots()
                fig.set_size_inches(8.0, 6.0)

                #
                df_no_duplicates = df.drop_duplicates(
                    subset=[column_1, column_2], keep="first",
                    ignore_index=True
                )

                #
                df_pivot = pd.pivot_table(
                    df_no_duplicates,
                    values=pipeline_component.criterion.name,
                    index=[column_1],
                    columns=[column_2],
                )
                df_pivot = df_pivot.sort_index(axis=0,
                                               ascending=False)
                sns.heatmap(
                    df_pivot,
                    ax=ax,
                    annot=True,
                    square=True,
                    fmt=".2f",
                    center=0,
                    cmap="turbo",
                    robust=True,
                    cbar_kws={"label": pipeline_component.criterion.name},
                )
                fig.autofmt_xdate(rotation=25)
                plt.yticks(rotation=0)

                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)
