import os

from pydock3.criterion.enrichment import __file__ as ENRICHMENT_MODULE_PATH


RANDOM_DATA_DIR_PATH = os.path.join(ENRICHMENT_MODULE_PATH, "random_classifier_probability")


def get_bonferroni_correction(
        n_actives: int,
        n_configurations: int,
        signif_level: float = 0.01,
        tables_dir: str = RANDOM_DATA_DIR_PATH,
):
    """
    :param n_actives: the number of actives in the retrospective dataset
    :param n_configurations: the number of docking configurations tested (i.e. number of different sets of parameters)
    :param signif_level: desired significance level of the test, .01 by default (who wants to be wrong 1/20 of the time)
    :param tables_dir: the directory where the tables are
    :return: the enrichment score threshold above which the whole endeavour has beaten random at p < signif_level
    """

    assert isinstance(n_actives, int)
    assert n_actives >= 1

    if n_actives > 50:
        print("Warning: no tables available for {} actives, reverting to 50 actives (which is more stringent).".format(
              n_actives))
        n_actives = 50

    threshold = float(signif_level / n_configurations)
    file_name = os.path.join(tables_dir, f"table_{n_actives}_actives.df")

    with open(file_name) as f:
        lines = f.readlines()[1:]

    found = False
    normalized_log_auc_thresh = None
    for line in lines:
        ll = line.split()
        pval = float(ll[-1])
        if pval <= threshold:
            found = True
            normalized_log_auc_thresh = float(ll[0])
            break

    if not found:
        raise ValueError("No threshold found (probably too few actives, or too many configurations)")

    return normalized_log_auc_thresh
