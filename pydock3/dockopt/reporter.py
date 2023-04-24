from __future__ import annotations
from typing import TYPE_CHECKING, NoReturn, Tuple, Union
import os
import itertools
import re
import base64

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import plotly.graph_objs as go

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

        #
        figures = []

        # Images
        for i, (_, row) in enumerate(pipeline_component.get_top_results_dataframe().head(1).iterrows()):
            best_job_dir_path = os.path.join(pipeline_component.best_retrodock_jobs_dir.path, f"rank={i+1}_step={row['component_id']}_conf={row['configuration_num']}")

            # ROC image
            roc_file_path = os.path.join(best_job_dir_path, ROC_PLOT_FILE_NAME)
            if os.path.exists(roc_file_path):
                figures.append(roc_file_path)

            # Energy plot image
            energy_plot_file_path = os.path.join(best_job_dir_path, ENERGY_PLOT_FILE_NAME)
            if os.path.exists(energy_plot_file_path):
                figures.append(energy_plot_file_path)

            # Energy plot image
            charge_plot_file_path = os.path.join(best_job_dir_path, CHARGE_PLOT_FILE_NAME)
            if os.path.exists(charge_plot_file_path):
                figures.append(charge_plot_file_path)

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
            figures.append(fig)

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
            figures.append(fig)

        # Generate the HTML report
        html = self.get_html(figures)
        with open(os.path.join(pipeline_component.component_dir.path, self.report_file_name), 'w') as f:
            f.write(html)

    @staticmethod
    def get_axis_label(param):
        label = ""
        barrier = ' ↳ '
        barrier_length = len(barrier)
        substrings = param.replace('_', ' ').split('.')
        max_len = max([len(s) + (barrier_length * i) for i, s in enumerate(substrings)])
        for i, s in enumerate(substrings):
            current_line = ""

            if i > 0:
                current_line += ' ' * barrier_length * (i - 1)
                current_line += barrier

            current_line += s
            current_line = current_line.ljust(max_len)
            current_line += '​'

            label += current_line
            if i < len(substrings) - 1:
                label += '<br>'

        return label

    @staticmethod
    def get_boxplot(df: pd.DataFrame, column: str, criterion_name: str, order=None, title=None) -> go.Figure:
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
                x=filtered_df[column].astype(str),  # treat as categorical
                name=str(value),
                boxpoints='all',  # Show all individual points
                jitter=0.6,  # Add some jitter for better visibility of points
                pointpos=-1.8,  # Position of points relative to the box
                whiskerwidth=0.5,
                marker=dict(color=c[i], size=6, line=dict(color='black', width=1)),  # Point marker color, size, and outline
                showlegend=False,
                offsetgroup=i,
            )

            # format the layout
            data.append(subplot)

        layout = go.Layout(
            font_family='monospace',
            xaxis=dict(
                title=HTMLReporter.get_axis_label(column.replace('parameters.', '')),
            ),
            yaxis=dict(
                title=HTMLReporter.get_axis_label(criterion_name),
            ),
            boxmode="overlay",
            boxgap=0.5,
        )
        fig = go.Figure(data=data, layout=layout)

        #
        fig.update_xaxes(
            type='category',
        )

        #
        if title is not None:
            fig.update_layout(title_text=title, title_xanchor="auto")

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
        heatmap = go.Heatmap(x=x_coords, y=y_coords, z=grid_scores, colorscale='Turbo', showscale=True)

        # Create the scatter plot
        scatter = go.Scatter(x=points[:, 0], y=points[:, 1], mode='markers',
                             marker=dict(color=scores_array, colorscale='Turbo', size=8,
                                         line=dict(color='gray', width=1)))

        # Combine the heatmap and scatter plot
        layout = go.Layout(
            font_family='monospace',
            xaxis=dict(title=HTMLReporter.get_axis_label(x.replace('parameters.', ''))),
            yaxis=dict(title=HTMLReporter.get_axis_label(y.replace('parameters.', ''))),
        )
        fig = go.Figure(data=[heatmap, scatter], layout=layout)

        #
        if title is not None:
            fig.update_layout(title_text=title, title_xanchor="auto")

        return fig

    @staticmethod
    def get_html(figures: Union[go.Figure, str]):
        divs = []

        for i, fig in enumerate(figures):
            is_one_slider = bool(isinstance(fig, str))
            if is_one_slider:
                image_path = fig
                with open(image_path, "rb") as image_file:
                    encoded_image = base64.b64encode(image_file.read()).decode("utf-8")

                slider_html = f"""
                    <label for="sizeSlider{i}">Size:</label>
                    <input type="range" min="100" max="1000" value="300" class="slider" id="sizeSlider{i}">
                    """
                divs.append(f"""<div>
                    {slider_html}
                    <img id="image{i}" src="data:image/jpeg;base64,{encoded_image}" style="max-width: 100%; max-height: 100%;">
                </div>
                """)
            else:
                plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
                slider_html = f"""
                    <label for="heightSlider{i}">Height:</label>
                    <input type="range" min="100" max="1000" value="300" class="slider" id="heightSlider{i}">
                    <label for="widthSlider{i}">Width:</label>
                    <input type="range" min="100" max="1000" value="500" class="slider" id="widthSlider{i}">
                    """

                divs.append(f"""<div>
                {slider_html}
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
            function setDefaultSize(imageElement, sizeSlider, widthSlider, heightSlider) {{
                if (sizeSlider) {{
                    var maxSize = Math.min(window.innerHeight, window.innerWidth) * 0.5;
                    sizeSlider.value = maxSize;
                    var aspectRatio = imageElement.naturalWidth / imageElement.naturalHeight;
                    imageElement.style.width = maxSize + 'px';
                    imageElement.style.height = maxSize / aspectRatio + 'px';
                }} else {{
                    var defaultHeight = window.innerHeight * 0.5;
                    var defaultWidth = window.innerWidth * 0.5;
                    heightSlider.value = defaultHeight;
                    widthSlider.value = defaultWidth;
                }}
            }}

            for (var i = 0; i < {len(figures)}; i++) {{
                var is_one_slider = document.getElementById('sizeSlider' + i) !== null;

                if (is_one_slider) {{
                var imageElement = document.getElementById('image' + i);
                var sizeSlider = document.getElementById('sizeSlider' + i);

                sizeSlider.oninput = (function(img) {{
                    return function() {{
                        var size = this.value;
                        var aspectRatio = img.naturalWidth / img.naturalHeight;
                        img.style.width = size + 'px';
                        img.style.height = size / aspectRatio + 'px';
                    }};
                }})(imageElement);

                setDefaultSize(imageElement, sizeSlider, null, null);
                }} else {{
                    var plotDiv = document.getElementById('plot' + i).getElementsByClassName('plotly-graph-div')[0];
                    var heightSlider = document.getElementById('heightSlider' + i);
                    var widthSlider = document.getElementById('widthSlider' + i);

                    heightSlider.oninput = (function(plot) {{
                        return function() {{
                            var height = this.value;
                            Plotly.update(plot, {{}}, {{height: height}});
                        }};
                    }})(plotDiv);
        
                    widthSlider.oninput = (function(plot) {{
                        return function() {{
                            var width = this.value;
                            Plotly.update(plot, {{}}, {{width: width}});
                        }};
                    }})(plotDiv);
        
                    setDefaultSize(null, null, widthSlider, heightSlider);
                }}
            }}
        </script>
    </body>
    </html>
        """
    
        return html
    
