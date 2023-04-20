from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn
import os
import itertools
import re
import base64
from io import BytesIO

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from PIL import Image
import plotly.graph_objs as go
import plotly.express as px

from pydock3.util import get_ordinal, sort_list_by_another_list
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

    # get the x dist and y dist per unit based on the two closest points
    sorted_coords = np.sort(np.unique(coords))
    min_distance = np.min(np.diff(sorted_coords))

    # Divide min_distance by min_units_between to get the increment value
    increment = min_distance / (min_units_between + 1)

    for i in range(len(sorted_coords) - 1):
        new_coords.extend(
            np.arange(sorted_coords[i], sorted_coords[i + 1], increment)
        )
    new_coords.append(sorted_coords[-1])

    return np.array(new_coords)


class HTMLReporter(Reporter):
    def __init__(self, report_file_name: str = "report.html") -> None:
        super().__init__(report_file_name)

    def write_report(self, pipeline_component: PipelineComponent) -> None:
        df = pipeline_component.load_results_dataframe()

        # Prepare data for subplots
        subplots = []
        subplot_titles = []

        # ROC Images
        for i, (_, row) in enumerate(pipeline_component.get_top_results_dataframe().iterrows()):
            best_job_dir_path = os.path.join(pipeline_component.best_retrodock_jobs_dir.path, f"{i + 1}_id={row['configuration_num']}")
            roc_file_path = os.path.join(best_job_dir_path, ROC_PLOT_FILE_NAME)
            if os.path.exists(roc_file_path):
                img = Image.open(roc_file_path)
                img_bytes = BytesIO()
                img.save(img_bytes, format='PNG')
                img_str = base64.b64encode(img_bytes.getvalue()).decode('utf-8')

                subplot_titles.append(f"linear-log ROC plot of {get_ordinal(i + 1)} best job<br>{best_job_dir_path}")
                subplots.append(go.Image(z=img_str))

        # Create boxplots
        boxplot_columns = sorted([
            column
            for column in df.columns
            if df[column].dropna().nunique() > 1  # multivalued parameters only, excluding NaNs
                and (
                    (column.startswith("parameters.dock_files_generation") or column.startswith("parameters.indock_file_generation"))  # include generation parameters
                    or column == "component_id"  # include component identifier
                )
        ])

        for column in boxplot_columns:
            fig = self.get_boxplot(df, column, pipeline_component.criterion.name)
            subplots.append(fig)
            subplot_titles.append(column.lstrip('parameters.').replace('.', '\n'))

        # Create heatmap for pairs of columns
        heatmap_numeric_columns = sorted([
            column
            for column in df.columns
            if df[column].dropna().nunique() > 1  # multivalued parameters only, excluding NaNs
                and (column.startswith("parameters.dock_files_generation") or column.startswith("parameters.indock_file_generation"))  # include generation parameters
                and pd.api.types.is_numeric_dtype(df[column])  # numeric parameters only
        ])
        df = df.sort_values(
            by=pipeline_component.criterion.name, ascending=False,
            ignore_index=True
        )

        for column_1, column_2 in itertools.combinations(heatmap_numeric_columns, 2):
            df_no_duplicates = df.drop_duplicates(subset=[column_1, column_2], keep="first", ignore_index=True)
            fig = self.get_heatmap(df_no_duplicates, column_1, column_2, pipeline_component.criterion.name)
            subplots.append(fig)
            subplot_titles.append(f"{column_2} vs. {column_1}")

        # Generate the HTML report
        html = self.get_html(subplots, subplot_titles)

        with open(os.path.join(pipeline_component.component_dir.path, self.report_file_name), 'w') as f:
            f.write(html)

    @staticmethod
    def get_boxplot(df: pd.DataFrame, column: str, criterion_name: str, order=None) -> go.Figure:
        if order is None:
            order = sorted(df[column].unique())

        # generate an array of rainbow colors by fixing the saturation and lightness of the HSL
        # representation of colour and marching around the hue.
        # Plotly accepts any CSS color format, see e.g. http://www.w3schools.com/cssref/css_colors_legal.asp.
        c = ['hsl(' + str(h) + ',50%' + ',50%)' for h in np.linspace(0, 360, len(order))]

        data = []
        for i, value in enumerate(order):
            filtered_df = df[df[column] == value]
            subplot = go.Box(
                y=filtered_df[criterion_name],
                x=filtered_df[column],
                name=str(value),
                boxpoints='all',  # Show all individual points
                jitter=0.5,  # Add some jitter for better visibility of points
                pointpos=-1.8,  # Position of points relative to the box
                whiskerwidth=0.3,
                marker=dict(color=c[i], size=6, line=dict(color='black', width=1)),  # Point marker color, size, and outline
                showlegend=False,
            )

            # format the layout
            data.append(subplot)

        layout = go.Layout(
            xaxis=dict(title=column),
            yaxis=dict(title=criterion_name),
            boxmode="group"
        )
        fig = go.Figure(data=data, layout=layout)
        return fig

    @staticmethod
    def get_heatmap(df: pd.DataFrame, x: str, y: str, scores: str, min_units_between: int = 20,
                    interp_method: str = 'cubic', title=None) -> go.Figure:
        # Extract points and scores
        points = df[[x, y]].to_numpy()
        scores_array = df[scores].to_numpy()

        # Create new x and y coordinates with the specified minimum number of grid units between any two points
        x_coords = create_new_coords(points[:, 0], min_units_between)
        y_coords = create_new_coords(points[:, 1], min_units_between)

        # Create mesh grid
        grid_x, grid_y = np.meshgrid(x_coords, y_coords)

        # Interpolate the data onto the grid
        grid_scores = griddata(points, scores_array, (grid_x, grid_y), method=interp_method)

        # Create the heatmap
        heatmap = go.Heatmap(x=x_coords, y=y_coords, z=grid_scores, colorscale='Turbo', showscale=False)

        # Create the scatter plot
        scatter = go.Scatter(x=points[:, 0], y=points[:, 1], mode='markers',
                             marker=dict(color=scores_array, colorscale='Turbo', size=8,
                                         line=dict(color='gray', width=1)))

        # Combine the heatmap and scatter plot
        fig = go.Figure(data=[heatmap, scatter])

        #
        fig.update_layout(xaxis_title=x, yaxis_title=y)
        if title is not None:
            fig.update_layout(title_text=title)

        return fig

    @staticmethod
    def get_html(subplots, subplot_titles):
        divs = []
        for i, (subplot, subplot_title) in enumerate(zip(subplots, subplot_titles)):
            fig = go.Figure(subplot)
            fig.update_layout(title=subplot_title, margin=dict(l=20, r=20, t=40, b=20))
            plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
            divs.append(f"""<div>
        <label for="heightSlider{i}">Height:</label>
        <input type="range" min="100" max="1000" value="300" class="slider" id="heightSlider{i}">
        <label for="widthSlider{i}">Width:</label>
        <input type="range" min="100" max="1000" value="500" class="slider" id="widthSlider{i}">
        <div id="plot{i}">{plot_html}</div>
    </div>
"""
                        )

        html = f"""
<!DOCTYPE html>
<html>
<head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        .slider {{
            width: 100%;
        }}
    </style>
</head>
<body>
    {''.join(divs)}
    <script>
        function getMaxHeight() {{
            return window.innerHeight || document.documentElement.clientHeight || document.body.clientHeight;
        }}

        function getMaxWidth() {{
            return window.innerWidth || document.documentElement.clientWidth || document.body.clientWidth;
        }}

        function setDefaultSize(plotDiv, heightSlider, widthSlider) {{
            var defaultHeight = getMaxHeight() * 0.5;
            var defaultWidth = getMaxWidth() * 0.5;
            heightSlider.value = defaultHeight;
            widthSlider.value = defaultWidth;
            Plotly.update(plotDiv, {{}}, {{height: defaultHeight, width: defaultWidth}});
        }}

        var plots = [];
        for (var i = 0; i < {len(subplots)}; i++) {{
            var plotDiv = document.getElementById('plot' + i).getElementsByClassName('plotly-graph-div')[0];
            plots.push(plotDiv);

            var heightSlider = document.getElementById('heightSlider' + i);
            heightSlider.max = getMaxHeight();
            heightSlider.oninput = (function(idx) {{
                return function() {{
                    var plotDiv = plots[idx];
                    var height = document.getElementById('heightSlider' + idx).value;
                    var width = document.getElementById('widthSlider' + idx).value;
                    Plotly.update(plotDiv, {{}}, {{height: height, width: width}});
                }};
            }})(i);

            var widthSlider = document.getElementById('widthSlider' + i);
            widthSlider.max = getMaxWidth();
            widthSlider.oninput = (function(idx) {{
                return function() {{
                    var plotDiv = plots[idx];
                    var height = document.getElementById('heightSlider' + idx).value;
                    var width = document.getElementById('widthSlider' + idx).value;
                    Plotly.update(plotDiv, {{}}, {{height: height, width: width}});
                }};
            }})(i);

            setDefaultSize(plotDiv, heightSlider, widthSlider);
        }}
    </script>
</body>
</html>
"""

        return html
