import os
import logging

import pandas as pd

from pydock3.criterion.enrichment import __file__ as ENRICHMENT_MODULE_INIT_PATH

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
ENRICHMENT_MODULE_PATH = os.path.dirname(ENRICHMENT_MODULE_INIT_PATH)
RANDOM_DATA_DIR_PATH = os.path.join(ENRICHMENT_MODULE_PATH, "random_classifier_probability")
MAX_TABLE_N_ACTIVES = 100


def get_random_classifier_performance_data(
        n_actives: int,
        tables_dir: str = RANDOM_DATA_DIR_PATH,
) -> pd.DataFrame:
    """
    Retrieves a DataFrame containing performance data for a random classifier.
    The DataFrame has the following columns:
    - normalized_log_auc: An observed normalized LogAUC value produced by a random classifier.
    - density: The density of the observed value under the null hypothesis.
    - prop: The proportion of the observed value under the null hypothesis.
    - cumul: The cumulative proportion of the observed value under the null hypothesis.
    - pval: The p-value, indicating the statistical significance of the observed value.

    :param n_actives: The number of positive results in the dataset.
    :param tables_dir: The directory where the data tables are stored.
    :return: A DataFrame containing performance data for a random classifier.
    """

    #
    if not isinstance(n_actives, int):
        raise TypeError(f"n_actives must be an integer, not {type(n_actives)}")
    if n_actives < 1:
        raise ValueError(f"n_actives must be >= 1, not {n_actives}")

    #
    if n_actives > MAX_TABLE_N_ACTIVES:
        n_actives_for_calculation = MAX_TABLE_N_ACTIVES
        logger.warning(f"No table available for {n_actives} actives, reverting to {n_actives_for_calculation} actives for calculation (which is more strict).")
    else:
        n_actives_for_calculation = n_actives

    #
    df = pd.read_csv(f"{tables_dir}/table_{n_actives_for_calculation}_actives.df", sep=" ")

    return df


def get_bonferroni_correction(
        n_actives: int,
        n_configurations: int,
        signif_level: float = 0.01,
        tables_dir: str = RANDOM_DATA_DIR_PATH,
) -> float:
    """
    :param n_actives: the number of actives in the retrospective dataset
    :param n_configurations: the number of docking configurations tested (i.e. number of different sets of parameters)
    :param signif_level: desired significance level of the test, .01 by default (who wants to be wrong 1/20 of the time)
    :param tables_dir: the directory where the tables are
    :return: the normalized LogAUC threshold above which the whole endeavor has beaten random at p < signif_level
    """

    #
    df_random = get_random_classifier_performance_data(n_actives, tables_dir)

    #
    threshold = float(signif_level / n_configurations)  # Bonferroni correction
    valid_thresholds = df_random[df_random['pval'] <= threshold]['normalized_log_auc']

    #
    if valid_thresholds.empty:
        raise ValueError("No threshold found (either too few actives or too many docking configurations tested)")

    #
    normalized_log_auc_thresh = valid_thresholds.iloc[0]

    return normalized_log_auc_thresh
