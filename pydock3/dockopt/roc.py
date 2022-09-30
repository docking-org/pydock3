from dataclasses import dataclass
import math

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt


@dataclass
class Point:
    x: float
    y: float


class ROC(object):
    def __init__(self, booleans, indices, alpha=None):
        #
        self.booleans = [bool(b) for b in booleans]
        self.indices = indices

        #
        self.num_actives = len([b for b in self.booleans if b])
        self.num_decoys = len([b for b in self.booleans if not b])

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
                    x_coord = float(num_decoys_witnessed_so_far / self.num_decoys)
                    y_coord = float(num_actives_witnessed_so_far / self.num_actives)
                    x_coords.append(x_coord)
                    y_coords.append(y_coord)
                    break  # last decoy has been seen, end of ROC, disregard remaining actives
                num_actives_witnessed_so_far += 1
                last_bool_was_decoy = False
            else:
                if not last_bool_was_decoy:
                    x_coord = float(num_decoys_witnessed_so_far/self.num_decoys)
                    y_coord = float(num_actives_witnessed_so_far/self.num_actives)
                    x_coords.append(x_coord)
                    y_coords.append(y_coord)
                num_decoys_witnessed_so_far += 1
                last_bool_was_decoy = True
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.points = [Point(x_coord, y_coord) for x_coord, y_coord in zip(self.x_coords, self.y_coords)]

        #
        x_coords_for_interpolation = [0.0] + self.x_coords
        y_coords_for_interpolation = [0.0] + self.y_coords
        self.f = lambda w: float(interpolate.interp1d(x_coords_for_interpolation, y_coords_for_interpolation, kind='previous')(w))

        #
        self.enrichment_score = self._get_enrichment_score()

    def _get_enrichment_score(self):
        return (self._get_literal_area_under_roc_curve_with_log_scaled_x_axis() - (1 - self.alpha)) / (-np.log(self.alpha) - (1 - self.alpha))

    def _get_literal_area_under_roc_curve_with_log_scaled_x_axis(self):
        #
        weights = []
        y_values_of_intervals = []
        last_x_value = self.alpha
        last_y_value = self.f(last_x_value)
        for i, (current_x_value, current_y_value) in enumerate(zip(self.x_coords, self.y_coords)):
            if current_x_value <= self.alpha:
                continue

            current_y_value = self.f(current_x_value)
            if current_y_value == last_y_value and (i+1) != len(self.y_coords):  # add for last y no matter what
                continue

            y_values_of_intervals.append(last_y_value)
            weights.append(np.log(float(current_x_value/last_x_value)))
            last_x_value = current_x_value
            last_y_value = current_y_value

        return np.dot(weights, y_values_of_intervals)



    def plot(self, save_path):
        # make plot of ROC curve of actives vs. decoys with log-scaled x-axis
        fig, ax = plt.subplots()
        fig.set_size_inches(8.0, 8.0)
        x_coords_for_plot = [self.alpha] + self.x_coords
        y_coords_for_plot = [self.f(self.alpha)] + self.y_coords
        ax.step(x_coords_for_plot, y_coords_for_plot, where='post', label=f"enrichment_score: {round(self.enrichment_score, 2)}")
        ax.legend()

        # draw curve of random classifier for reference
        ax.semilogx([float(i / 1000) for i in range(0, 1001)], [float(i / 1000) for i in range(0, 1001)], "--", linewidth=1, color='Black')

        # set axis labels
        ax.set_xlabel('top fraction of decoys')
        ax.set_ylabel('top fraction of actives')

        # set log scale x-axis
        ax.set_xscale('log')

        # set plot axis limits
        ax.set_xlim(left=self.alpha, right=1.0)
        ax.set_ylim(bottom=0.0, top=1.0)

        # set axis ticks
        order_of_magnitude = math.floor(math.log(self.num_decoys, 10))
        ax.set_xticks([self.alpha] + [float(10 ** x) for x in range(-order_of_magnitude, 1, 1)])
        ax.set_yticks([float(j / 10) for j in range(0, 11)])

        # set title and include alpha to account for inconsistent inclusion of alpha in xticks
        ax.set_title(f"Log ROC (num_decoys={self.num_decoys}, num_actives={self.num_actives}, cutoff={np.format_float_scientific(self.alpha, precision=3)})")

        # save image and close
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close(fig)
