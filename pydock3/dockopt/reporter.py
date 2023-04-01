from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os
import itertools
import re

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


from pydock3.util import get_ordinal, sort_list_by_another_list
from pydock3.files import File
from pydock3.retrodock.retrodock import ROC_PLOT_FILE_NAME, ENERGY_PLOT_FILE_NAME, CHARGE_PLOT_FILE_NAME
if TYPE_CHECKING:
    from pydock3.dockopt.pipeline import PipelineComponent


class Reporter(object):
    def __init__(self, report_file_name: str) -> None:
        self.report_file_name = report_file_name

    def write_report(self, pipeline_component: PipelineComponent) -> NoReturn:
        raise NotImplementedError


def get_sorted_component_ids(component_ids: list) -> list:
    """Returns a list of component ids sorted by the order in which they should be processed."""

    numeric_ids = [
        tuple([
            int(re.sub('\D', '', piece))  # remove non-numeric characters, cast to int
            for piece in component_id.split('.')  # split on periods (e.g. '1.2.3' -> ['1', '2', '3']
        ])
        for component_id in component_ids
    ]

    return sort_list_by_another_list(component_ids, numeric_ids)  # sort by numeric ids, return component ids


def create_new_coords(coords: np.ndarray, min_units_between: int) -> np.ndarray:
    new_coords = []

    sorted_coords = np.sort(np.unique(coords))

    for i in range(len(sorted_coords) - 1):
        new_coords.extend(
            np.linspace(sorted_coords[i], sorted_coords[i + 1], num=min_units_between + 1, endpoint=False))
    new_coords.append(sorted_coords[-1])

    return np.array(new_coords)


def heatmap(df: pd.DataFrame, x: str, y: str, scores: str, min_units_between: int = 5, interp_method: str = 'nearest') -> None:
    # Extract points and scores

    points = df[[x, y]].to_numpy()
    scores_array = df[scores].to_numpy()

    # Create new x and y coordinates with the specified minimum number of grid units between points
    x_coords = create_new_coords(points[:, 0], min_units_between)
    y_coords = create_new_coords(points[:, 1], min_units_between)

    # Create the uneven mesh grid
    grid_x, grid_y = np.meshgrid(x_coords, y_coords)

    # Interpolate the data onto the grid
    grid_scores = griddata(points, scores_array, (grid_x, grid_y), method=interp_method, fill_value=0)

    # Define colormap and normalization for both heatmap and scatter plot
    cmap = plt.get_cmap('hot')
    norm = Normalize(vmin=scores_array.min(), vmax=scores_array.max())

    # Plot the heatmap
    heatmap_image = plt.imshow(grid_scores, origin='lower',
                               extent=(x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()),
                               cmap=cmap, aspect='auto', norm=norm)
    cbar = plt.colorbar(heatmap_image, label='Scores')

    # Plot the scatter points with the same colormap and normalization
    plt.scatter(points[:, 0], points[:, 1], c=scores_array, cmap=cmap, edgecolors='k', norm=norm)

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title('Heatmap of Scores')
    plt.show()


class PDFReporter(Reporter):
    def __init__(self, report_file_name: str = "report.pdf") -> None:
        super().__init__(report_file_name)

        #
        plt.rcParams.update({"font.size": 14})

    def write_report(self, pipeline_component: PipelineComponent) -> None:

        df = pipeline_component.load_results_dataframe()

        with PdfPages(os.path.join(pipeline_component.component_dir.path, self.report_file_name)) as f:
            #
            fig = plt.figure(figsize=(11.0, 8.5))

            #
            for i, (_, row) in enumerate(pipeline_component.get_top_results_dataframe().iterrows()):
                #
                best_job_dir_path = os.path.join(pipeline_component.best_retrodock_jobs_dir.path, f"{i + 1}_id={row['configuration_num']}")

                #
                roc_file_path = os.path.join(
                    best_job_dir_path, ROC_PLOT_FILE_NAME
                )
                if File.file_exists(roc_file_path):  # TODO: replace this logic with call to `Retrodock.analyze` method or something in order to ensure that these Retrodock-pertaining files are created
                    image = mpimg.imread(roc_file_path)
                    plt.axis("off")
                    plt.suptitle(f"linear-log ROC plot of {get_ordinal(i + 1)} best job\n{best_job_dir_path}")
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
                    plt.suptitle("energy terms ridgeline plot of best job")
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
            boxplot_columns = [
                column
                for column in df.columns
                if df[column].nunique() > 1  # multivalued parameters only
                    and (
                        (column.startswith("parameters.dock_files_generation") or column.startswith("parameters.indock_file_generation"))  # include generation parameters
                        or column == "component_id"  # include component identifier
                    )
            ]
            for column in boxplot_columns:
                fig = plt.figure(figsize=(8, 6))

                # need to sort component ids by order of numerical pieces
                if column == "component_id":
                    order = get_sorted_component_ids(list(df["component_id"].unique()))
                else:
                    order = None  # seaborn will infer order from data

                #
                if len(df) == len(df.drop_duplicates(subset=[column])):  # use bar plot instead if all values are unique
                    sns.barplot(data=df, x=column, y=pipeline_component.criterion.name, order=order, rot=0)
                else:
                    sns.boxplot(
                        data=df,
                        x=column,
                        y=pipeline_component.criterion.name,
                        order=order,
                        showfliers=False,
                        boxprops={"facecolor": "None"},
                    )
                    sns.stripplot(data=df, x=column, y=pipeline_component.criterion.name, order=order, zorder=0.5)

                #
                fig.autofmt_xdate(rotation=25)
                plt.yticks(rotation=0)
                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)

            #
            heatmap_columns = [
                column
                for column in df.columns
                if df[column].nunique() > 1  # multivalued parameters only
                    and (column.startswith("parameters.dock_files_generation") or column.startswith("parameters.indock_file_generation"))  # include generation parameters
                    and pd.api.types.is_numeric_dtype(df[column])  # numeric parameters only
            ]
            df = df.sort_values(
                by=pipeline_component.criterion.name, ascending=False,
                ignore_index=True
            )
            for column_1, column_2 in itertools.combinations(
                    heatmap_columns, 2
            ):
                #
                fig, ax = plt.subplots()
                fig.set_size_inches(8.0, 6.0)

                #
                df_no_duplicates = df.drop_duplicates(subset=[column_1, column_2], keep="first", ignore_index=True)

                #
                """
                df_pivot = pd.pivot_table(
                    df_no_duplicates,
                    values=pipeline_component.criterion.name,
                    index=[column_1],
                    columns=[column_2],
                )
                df_pivot = df_pivot.sort_index(axis=0, ascending=False)
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
                """

                #
                heatmap(df_no_duplicates, column_1, column_2, pipeline_component.criterion.name, min_units_between=5, interp_method='nearest')

                #
                fig.autofmt_xdate(rotation=25)
                plt.yticks(rotation=0)
                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)
