from collections import namedtuple
import sys

import numpy as np
from scipy import integrate, interpolate, special
from pytoc import TOC


#
DEFAULT_ALPHA = 1.e-3


#
class EnrichmentMetric(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    @staticmethod
    def calculate(*args, **kwargs):
        raise NotImplementedError


#
class LogAUC(EnrichmentMetric):

    def __init__(self, name, f, alpha, discontinuous_x_coords=None):
        self.f = f
        self.alpha = alpha
        self.discontinuous_x_coords = discontinuous_x_coords
        value = self.calculate(self.f, self.alpha, discontinuous_x_coords=self.discontinuous_x_coords)

        super().__init__(name, value)

    @staticmethod
    def calculate(f, alpha, discontinuous_x_coords=None):
        return LogAUC.raw_log_auc(f, alpha, discontinuous_x_coords=discontinuous_x_coords) / LogAUC.raw_log_auc(best_roc_curve_function, alpha)

    @staticmethod
    def raw_log_auc(f, alpha, subintervals_limit=int(10e6), discontinuous_x_coords=None, abs_err_tol=1e-4):
        """This is the literal 'area under the curve' of the given function f over a logarithmic scale on the x-axis."""
        
        min_alpha = 1.e-323  # any smaller value results in divide by zero encountered in np.log(alpha)
        max_alpha = 1.0  # alpha >= 1.0 breaks integral
        if not (min_alpha <= alpha < max_alpha):
            raise ValueError(f"Violates constraint: {min_alpha} <= alpha < {max_alpha}")
        return integrate.quad(lambda x: f(np.exp(x)), np.log(alpha), 0, limit=subintervals_limit,
                              points=discontinuous_x_coords, epsabs=abs_err_tol)[0]


class EnrichmentScore(EnrichmentMetric):

    def __init__(self, name, f, alpha, discontinuous_x_coords=None):
        self.f = f
        self.alpha = alpha
        self.discontinuous_x_coords = discontinuous_x_coords
        value = self.calculate(self.f, self.alpha, discontinuous_x_coords=self.discontinuous_x_coords)

        super().__init__(name, value)

    @staticmethod
    def calculate(f, alpha, discontinuous_x_coords=None):
        best_roc_curve_auc_normalized = LogAUC.calculate(best_roc_curve_function, alpha)
        random_roc_curve_auc_normalized = LogAUC.calculate(random_roc_curve_function, alpha)
        return (LogAUC.calculate(f, alpha, discontinuous_x_coords=discontinuous_x_coords) - random_roc_curve_auc_normalized) / (best_roc_curve_auc_normalized - random_roc_curve_auc_normalized)


#
class UnbalancedLogAUC(LogAUC):
    METRIC_NAME = "unbalanced_log_auc"

    def __init__(self, f, discontinuous_x_coords=None):

        super().__init__(name=self.METRIC_NAME, f=f, alpha=DEFAULT_ALPHA, discontinuous_x_coords=discontinuous_x_coords)


class BalancedLogAUC(LogAUC):
    METRIC_NAME = "balanced_log_auc"

    def __init__(self, f, num_x, discontinuous_x_coords=None):
        alpha = self.balanced_alpha(num_x=num_x)
        super().__init__(name=self.METRIC_NAME, f=f, alpha=alpha, discontinuous_x_coords=discontinuous_x_coords)

    @staticmethod
    def balanced_alpha(num_x):
        return 1.0 / (2.0 * num_x)


class UnbalancedEnrichmentScore(EnrichmentScore):
    METRIC_NAME = "unbalanced_enrichment_score"

    def __init__(self, f, discontinuous_x_coords=None):
        super().__init__(name=self.METRIC_NAME, f=f, alpha=DEFAULT_ALPHA, discontinuous_x_coords=discontinuous_x_coords)


class BalancedEnrichmentScore(EnrichmentScore):
    METRIC_NAME = "balanced_enrichment_score"

    def __init__(self, f, num_x, discontinuous_x_coords=None):
        alpha = self.balanced_alpha(num_x=num_x)
        super().__init__(name=self.METRIC_NAME, f=f, alpha=alpha, discontinuous_x_coords=discontinuous_x_coords)

    @staticmethod
    def balanced_alpha(num_x, log_base=np.e):
        return -special.lambertw(-np.log(log_base) / (2.0 * num_x)).real / np.log(log_base)


#
def best_roc_curve_function(x):
    return 1


def worst_roc_curve_function(x):
    return 0


def random_roc_curve_function(x):
    return x


#
def is_sorted(l):
    return all(l[i] <= l[i + 1] for i in range(len(l) - 1))


def step_function(x, y):
    #
    if any([not (0.0 <= i <= 1.0) for i in x]):
        raise ValueError("Violates constraint: 0.0 <= x <= 1.0")
    if any([not (0.0 <= i <= 1.0) for i in y]):
        raise ValueError("Violates constraint: 0.0 <= y <= 1.0")

    #
    x, y = zip(*sorted(zip(x, y)))  # sort both x and y by x
    if not (is_sorted(x) and is_sorted(y)):  # not monotonically increasing
        raise ValueError("x and y must be taken from a monotonically increasing function.")
    
    #
    return lambda w: float(interpolate.interp1d(x, y, kind='previous')(w))


#
EnrichmentAnalysis = namedtuple("EnrichmentAnalysis", "unbalanced_log_auc balanced_log_auc unbalanced_enrichment_score balanced_enrichment_score")


def get_enrichment_analysis(x, y):
    f = step_function(x, y)
    discontinuous_x_coords = list(set([0.0] + x + [1.0]))

    return EnrichmentAnalysis(
        unbalanced_log_auc=UnbalancedLogAUC(f=f, discontinuous_x_coords=discontinuous_x_coords),
        balanced_log_auc=BalancedLogAUC(f=f, discontinuous_x_coords=discontinuous_x_coords, num_x=len(x)-1),  # num interdecoy intervals
        unbalanced_enrichment_score=UnbalancedEnrichmentScore(f=f, discontinuous_x_coords=discontinuous_x_coords),
        balanced_enrichment_score=BalancedEnrichmentScore(f=f, discontinuous_x_coords=discontinuous_x_coords, num_x=len(x)-1),  # num interdecoy intervals
    )


#
def get_toc_points(booleans, indices):
    booleans_array = np.array(booleans)
    indices_array = np.array(indices)
    thresholds_array = np.unique(indices_array)
    toc = TOC(booleans_array, indices_array, thresholds_array)
    toc_points_x = np.array(toc.TOCX[0], dtype=float)
    toc_points_y = np.array(toc.TOCY[0], dtype=float)
    toc_points_x, toc_points_y = zip(*sorted(zip(toc_points_x, toc_points_y)))
    toc_points_x = np.array(toc_points_x)
    toc_points_y = np.array(toc_points_y)

    return toc_points_x, toc_points_y


#
def get_roc_points(booleans, indices):
    # get TOC points
    toc_points_x, toc_points_y = get_toc_points(booleans, indices)

    # get number of items in each binary class
    num_positive = len([b for b in booleans if b])
    num_negative = len([b for b in booleans if not b])

    # get ROC points
    roc_points_x = toc_points_x - toc_points_y  # unshear the ROC
    roc_points_y = toc_points_y
    mask = np.zeros(len(roc_points_x), dtype=bool)
    roc_points_x_unique = list(np.unique(roc_points_x))
    for x_unique in roc_points_x_unique:
        indices = [index for index, x in enumerate(roc_points_x) if x == x_unique]
        # set mask element to True only if index corresponds to max y for given x
        max_index = indices[0]
        for index in indices:
            if toc_points_y[index] > toc_points_y[max_index]:
                max_index = index
        mask[max_index] = True
    roc_points_x = roc_points_x[mask] / num_negative
    roc_points_y = roc_points_y[mask] / num_positive

    return roc_points_x, roc_points_y
