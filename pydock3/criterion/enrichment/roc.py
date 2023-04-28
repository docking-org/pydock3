from typing import Iterable
from dataclasses import dataclass
import math

import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')  # set the backend to Agg (no interactive plots)
from matplotlib import pyplot as plt

plt.rcParams.update({"font.size": 14})


@dataclass
class Point:
    x: float
    y: float


class ROC(object):
    def __init__(
        self,
        booleans: Iterable[bool],
        alpha: float = None,
    ):
        #
        self.booleans = booleans

        #
        self.num_actives = len([b for b in self.booleans if b])
        self.num_decoys = len([b for b in self.booleans if not b])

        # validate num actives and num decoys
        if self.num_actives == 0 or self.num_decoys == 0:
            raise ValueError(f"Number of actives and number of decoys both must be greater than zero!\n\tnum_actives={self.num_actives}\n\tnum_decoys={self.num_decoys}")

        # validate and set alpha
        if alpha is None:
            alpha = float(1 / (self.num_decoys * np.e))
        if not ((alpha > 0.0) and (alpha < 1.0)):
            raise ValueError("ROC alpha must be in range (0, 1)")
        self.alpha = alpha

        # get ROC points
        x_coords = []
        y_coords = []
        num_decoys_witnessed_so_far = 0
        num_actives_witnessed_so_far = 0
        last_bool_was_decoy = False
        for b in self.booleans:
            if b:
                if num_decoys_witnessed_so_far == self.num_decoys:
                    break  # last decoy has been seen, end of ROC, disregard remaining actives
                num_actives_witnessed_so_far += 1
                last_bool_was_decoy = False
            else:
                if not last_bool_was_decoy:
                    x_coord = float(num_decoys_witnessed_so_far / self.num_decoys)
                    y_coord = float(num_actives_witnessed_so_far / self.num_actives)
                    x_coords.append(x_coord)
                    y_coords.append(y_coord)
                num_decoys_witnessed_so_far += 1
                last_bool_was_decoy = True
        self.x_coords = x_coords  # num points = num_decoys, each ith point represents interval [i/n, (i+1)/n]
        self.y_coords = y_coords
        self.points = [
            Point(x_coord, y_coord)
            for x_coord, y_coord in zip(self.x_coords, self.y_coords)
        ]

        #
        x_coords_for_interpolation = self.x_coords + [1.0]  # add point at x=1.0 to complete the last interval [(n-1)/n, 1.0]
        y_coords_for_interpolation = self.y_coords + [self.y_coords[-1]]
        self.f = lambda w: float(
            interpolate.interp1d(
                x_coords_for_interpolation, y_coords_for_interpolation, kind="previous"
            )(w)
        )

        #
        self._literal_log_auc = self._get_literal_log_auc()
        self._random_literal_log_auc = self._get_random_literal_log_auc()
        self._optimal_literal_log_auc = self._get_optimal_literal_log_auc()
        #self.log_auc = self._get_log_auc()  # unnormalized LogAUC should probably be avoided entirely
        self.normalized_log_auc = self._get_normalized_log_auc()

    def _get_random_literal_log_auc(self) -> float:
        return float(1 - self.alpha)
    
    def _get_optimal_literal_log_auc(self) -> float:
        return -np.log(self.alpha)

    def _get_literal_log_auc(self) -> float:
        # remove point at x=0.0 and add point at x=1.0
        x_values = self.x_coords[1:] + [1.0]
        y_values = self.y_coords[1:] + [self.y_coords[-1]]

        #
        weights = []
        y_values_of_intervals = []
        previous_x_value = self.alpha
        previous_y_value = self.f(self.alpha)  # initialize according to point at x=alpha
        for i, (current_x_value, current_y_value) in enumerate(
            zip(x_values, y_values)
        ):
            if current_y_value == previous_y_value and (i + 1) != len(
                y_values
            ):  # add for last y no matter what
                continue

            #
            y_values_of_intervals.append(previous_y_value)
            weights.append(np.log(float(current_x_value / previous_x_value)))

            #
            previous_x_value = current_x_value
            previous_y_value = current_y_value

        return float(np.dot(weights, y_values_of_intervals).item())
    
    def _get_log_auc(self) -> float:
        return self._literal_log_auc / self._optimal_literal_log_auc
    
    def _get_normalized_log_auc(self) -> float:
        return (self._literal_log_auc - self._random_literal_log_auc) / (self._optimal_literal_log_auc - self._random_literal_log_auc)

    def plot(
        self,
        save_path: str
    ) -> None:
        fig, ax = plt.subplots()
        fig.set_size_inches(8.0, 8.0)

        # draw curve of random classifier for reference
        ax.semilogx(
            [float(i / 1000) for i in range(0, 1001)],
            [float(i / 1000) for i in range(0, 1001)],
            "--",
            linewidth=1,
            color="Black",
            label=f"random classifier",
        )

        # make plot of ROC curve of actives vs. decoys with log-scaled x-axis
        x_coords_for_plot = self.x_coords + [1.0]  # add point at x=1.0 to complete the last interval [(n-1)/n, 1.0]
        y_coords_for_plot = self.y_coords + [self.y_coords[-1]]
        ax.step(
            x_coords_for_plot,
            y_coords_for_plot,
            where="post",
            label=f"ROC curve",
        )

        # add an extra label of normalized LogAUC (nothing extra will be plotted)
        plt.plot([], [], " ", label=f"# of actives: {self.num_actives}")
        plt.plot([], [], " ", label=f"# of decoys: {self.num_decoys}")
        plt.plot(
            [],
            [],
            " ",
            label=f"x-axis interval: [{np.format_float_scientific(self.alpha, precision=3)}, 1.0]",
        )
        plt.plot(
            [], [], " ", label=f"normalized LogAUC: {round(self.normalized_log_auc, 3)}"
        )

        # set legend
        ax.legend(framealpha=0.75)

        # set axis labels
        ax.set_xlabel("false positive rate (i.e., top fraction of decoys accepted)")
        ax.set_ylabel("true positive rate (i.e., top fraction of actives accepted)")

        # set log scale x-axis
        ax.set_xscale("log")

        # set plot axis limits
        ax.set_xlim(left=self.alpha, right=1.0)
        ax.set_ylim(bottom=0.0, top=1.0)

        # set axis ticks
        order_of_magnitude = -math.floor(math.log(self.alpha, 10)) - 1
        ax.set_xticks(
            [self.alpha] + [float(10**x) for x in range(-order_of_magnitude, 1, 1)]
        )
        ax.set_yticks([float(j / 10) for j in range(0, 11)])

        # set title
        ax.set_title(f"linear-log ROC plot")

        # save image and close
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig)
