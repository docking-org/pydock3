import os
import logging
from uuid import uuid4
import time
from dataclasses import astuple
from typing import List, Union, Tuple
import collections

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joypy import joyplot

from pydock3.util import Script
from pydock3.files import (
    Dir,
    File,
    IndockFile,
    OutdockFile,
)
from pydock3.retrodock.retrospective_dataset import RetrospectiveDataset
from pydock3.criterion.enrichment.roc import ROC
from pydock3.jobs import ArrayDockingJob, OUTDOCK_FILE_NAME
from pydock3.blastermaster.blastermaster import BlasterFiles, BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT
from pydock3.jobs import JobSubmissionResult
from pydock3.job_schedulers import SGEJobScheduler, SlurmJobScheduler
from pydock3.docking import __file__ as DOCKING_INIT_FILE_PATH

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
DOCK3_EXECUTABLE_PATH = os.path.join(
    os.path.dirname(DOCKING_INIT_FILE_PATH), "dock3", "dock64"
)

#
SCHEDULER_NAME_TO_CLASS_DICT = {
    "sge": SGEJobScheduler,
    "slurm": SlurmJobScheduler,
}

#
NORMALIZED_LOG_AUC_FILE_NAME = "normalized_log_auc.txt"
ROC_PLOT_FILE_NAME = "roc.png"
ENERGY_TERMS_PLOT_FILE_NAME = "energy.png"
CHARGE_PLOT_FILE_NAME = "charge.png"


BINARY_CLASS_COLOR_PALETTE = collections.OrderedDict({
    "positive": "#686de0",
    "negative": "#eb4d4b",
})


# TODO: add this as a decorator of `submit` method of `DockoptJob`
def log_job_submission_result(job, submission_result, procs):
    if submission_result is JobSubmissionResult.SUCCESS:
        logger.info(f"Job '{job.name}' successfully submitted.")
    elif submission_result is JobSubmissionResult.FAILED:
        for proc in procs:
            raise Exception(
                f"Job submission failed for '{job.name}' due to error: {proc.stderr}"
            )
    elif submission_result is JobSubmissionResult.SKIPPED_BECAUSE_ALREADY_COMPLETE:
        logger.info(
            f"Job submission skipped for '{job.name}' since all its OUTDOCK files already exist."
        )
    elif submission_result is JobSubmissionResult.SKIPPED_BECAUSE_STILL_ON_JOB_SCHEDULER_QUEUE:
        logger.info(
            f"Job submission skipped for '{job.name}' since it is still running from a previous submission."
        )
    else:
        raise Exception(f"Unrecognized JobSubmissionResult: {submission_result}")


def str_to_float(s, alternative_if_uncastable=np.nan):
    """cast numerical fields as float"""

    try:
        result = float(s)
    except ValueError:
        result = alternative_if_uncastable
    return result


def get_results_dataframe_from_positives_job_and_negatives_job_outdock_files(
    positives_outdock_file_path, negatives_outdock_file_path
):
    """build dataframe of docking results from outdock files"""

    #
    positives_outdock_file = OutdockFile(positives_outdock_file_path)
    negatives_outdock_file = OutdockFile(negatives_outdock_file_path)

    #
    positives_outdock_df = positives_outdock_file.get_dataframe()
    negatives_outdock_df = negatives_outdock_file.get_dataframe()

    # set is_positive column based on outdock file
    positives_outdock_df["is_positive"] = [1 for _ in range(len(positives_outdock_df))]
    negatives_outdock_df["is_positive"] = [0 for _ in range(len(negatives_outdock_df))]

    # set class_label column based on outdock file
    positives_outdock_df["class_label"] = [
        "positive" for _ in range(len(positives_outdock_df))
    ]
    negatives_outdock_df["class_label"] = [
        "negative" for _ in range(len(negatives_outdock_df))
    ]

    # build dataframe of docking results from outdock files
    df = pd.DataFrame()
    df = pd.concat([df, positives_outdock_df], ignore_index=True)
    df = pd.concat([df, negatives_outdock_df], ignore_index=True)

    # replace relevant str columns with float equivalents & change column names
    for old_col, new_col in [
        ("Total", "total_energy"),
        ("elect", "electrostatic_energy"),
        ("vdW", "vdw_energy"),
        ("psol", "polar_desolvation_energy"),
        ("asol", "apolar_desolvation_energy"),
        ("charge", "charge"),
    ]:
        df[new_col] = df[old_col].apply(lambda s: str_to_float(s))
        if new_col != old_col:
            df = df.drop(old_col, axis=1)

    # sort dataframe by total energy score
    df = df.sort_values(
        by=["total_energy", "is_positive"],
        na_position="last",
        ignore_index=True,
    )  # sorting secondarily by 'is_positive' (0 or 1) ensures that negatives are ranked before positives in case they have the same exact score (pessimistic approach)
    df = df.drop_duplicates(
        subset=["db2_file_path"],
        keep="first",
        ignore_index=True
    )  # keep only the best score per molecule

    return df


def make_ridgeline_plot_of_energy_terms(
        df: pd.DataFrame,
        save_path: Union[None, str] = None,
        title: Union[None, str] = None,
        figsize: Tuple[int, int] = (8, 8),
        dpi: int = 300,
        alpha: float = 0.8,
) -> Tuple[plt.Figure, plt.Axes]:
    """make ridgeline plot of energy terms"""

    #
    if save_path is None:
        save_path = ENERGY_TERMS_PLOT_FILE_NAME

    #
    if title is None:
        title = "Distributions of Energy Terms (Positives vs. Negatives)"

    #
    columns = [
        "total_energy",
        "electrostatic_energy",
        "vdw_energy",
        "polar_desolvation_energy",
        "apolar_desolvation_energy",
    ]

    #
    pivot_rows = []
    for i, row in df.iterrows():
        for col in columns:
            pivot_row = {"energy_term": col}
            if row["is_positive"] == 1:
                pivot_row["positive"] = str_to_float(row[col])
                pivot_row["negative"] = np.nan
            else:
                pivot_row["positive"] = np.nan
                pivot_row["negative"] = str_to_float(row[col])
            pivot_rows.append(pivot_row)
    df_pivot = pd.DataFrame(pivot_rows)
    fig, ax = joyplot(
        data=df_pivot,
        by="energy_term",
        column=["positive", "negative"],
        color=[BINARY_CLASS_COLOR_PALETTE['positive'], BINARY_CLASS_COLOR_PALETTE['negative']],
        legend=True,
        alpha=alpha,
        figsize=figsize,
        ylim="own",
    )
    ax[-1].set_xlabel("Δ Energy (kcal/mol)")  # set x-axis label on last subplot to ensure it appears at bottom of plot
    plt.title(title)
    try:
        plt.tight_layout()
    except UserWarning:
        pass
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    return fig, ax


def make_split_violin_plot_of_charge(
        df: pd.DataFrame,
        save_path: Union[None, str] = None,
        title: Union[None, str] = None,
        figsize: Tuple[int, int] = (8, 8),
        dpi: int = 300,
) -> Tuple[plt.Figure, plt.Axes]:
    """make split violin plot of charge"""

    #
    if save_path is None:
        save_path = CHARGE_PLOT_FILE_NAME

    #
    if title is None:
        title = 'Distributions of Charge (Positives vs. Negatives)'

    #
    fig, ax = plt.subplots(figsize=figsize)
    c = sns.color_palette(
        palette=[BINARY_CLASS_COLOR_PALETTE['positive'], BINARY_CLASS_COLOR_PALETTE['negative']],
        n_colors=2,
    )
    palette = {'positive': c[0], 'negative': c[1]}
    sns.violinplot(
        data=df,
        x="charge",
        y="total_energy",
        split=True,
        hue="class_label",
        inner="quartile",
        palette=palette,
    )
    plt.title(title)
    try:
        plt.tight_layout()
    except UserWarning:
        pass
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    return fig, ax


def process_retrodock_job_results(
        positives_outdock_file_path: str,
        negatives_outdock_file_path: str,
        save_dir_path: str,
):
    """process retrodock job results"""

    # set save file paths
    normalized_log_auc_save_path = os.path.join(save_dir_path, NORMALIZED_LOG_AUC_FILE_NAME)
    roc_plot_save_path = os.path.join(save_dir_path, ROC_PLOT_FILE_NAME)
    energy_terms_plot_save_path = os.path.join(save_dir_path, ENERGY_TERMS_PLOT_FILE_NAME)
    charge_plot_save_path = os.path.join(save_dir_path, CHARGE_PLOT_FILE_NAME)

    # load results
    df = get_results_dataframe_from_positives_job_and_negatives_job_outdock_files(
        positives_outdock_file_path,
        negatives_outdock_file_path,
    )

    # calculate ROC
    roc = ROC(booleans=df["is_positive"].astype(bool))
    with open(normalized_log_auc_save_path, "w") as f:
        f.write(f"{roc.normalized_log_auc}")

    # make plots
    roc.plot(save_path=roc_plot_save_path)
    make_ridgeline_plot_of_energy_terms(df, save_path=energy_terms_plot_save_path)
    make_split_violin_plot_of_charge(df, save_path=charge_plot_save_path)


class Retrodock(Script):
    JOB_DIR_NAME = "retrodock_job"
    DOCK_FILES_DIR_NAME = "dockfiles"
    INDOCK_FILE_NAME = "INDOCK"
    POSITIVES_TGZ_FILE_NAME = "positives.tgz"
    NEGATIVES_TGZ_FILE_NAME = "negatives.tgz"

    SINGLE_TASK_NUM = 1

    def __init__(self):
        super().__init__()

    def new(
            self,
            job_dir_path=JOB_DIR_NAME,
            overwrite=False  # TODO: see blastermaster script for example
    ) -> None:
        """create new job directory with dockfiles & INDOCK file"""

        # create job dir
        job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create dockfiles dir & copy in dock files & INDOCK
        dock_files_dir_path = os.path.join(job_dir.path, self.DOCK_FILES_DIR_NAME)
        dock_files_dir = Dir(dock_files_dir_path, create=True)
        dock_files = BlasterFiles(dock_files_dir).dock_files

        #
        dock_files_required = [dock_file for dock_file in astuple(dock_files) if dock_file.name != BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT['electrostatics_phi_size_file']]
        docking_configuration_file_names_required = [os.path.join(dock_files_dir.name, dock_file.name) for dock_file in dock_files_required] + [self.INDOCK_FILE_NAME]
        docking_configuration_file_names_not_detected = [file_name for file_name in docking_configuration_file_names_required if not File.file_exists(file_name)]
        if len(docking_configuration_file_names_not_detected) > 0:
            raise Exception("Missing required input files in current working directory:\n"+"\n\t".join(docking_configuration_file_names_not_detected))
        else:
            logger.info(
                f"Copying the following files from current working directory into job dockfiles directory:\n\t{docking_configuration_file_names_required}"
            )
            for file_name in docking_configuration_file_names_required:
                dock_files_dir.copy_in_file(file_name)

        # copy in positives and negatives TGZ files
        tgz_files = [self.POSITIVES_TGZ_FILE_NAME, self.NEGATIVES_TGZ_FILE_NAME]
        tgz_file_names_in_cwd = [f for f in tgz_files if os.path.isfile(f)]
        tgz_file_names_not_in_cwd = [f for f in tgz_files if not os.path.isfile(f)]
        if tgz_file_names_in_cwd:
            files_to_copy_str = "\n\t".join(tgz_file_names_in_cwd)
            logger.info(
                f"Copying the following files from current directory into job directory:\n\t{files_to_copy_str}"
            )
            for tgz_file_name in tgz_file_names_in_cwd:
                job_dir.copy_in_file(tgz_file_name)
        if tgz_file_names_not_in_cwd:
            files_missing_str = "\n\t".join(tgz_file_names_not_in_cwd)
            logger.info(
                f"The following required files were not found in current working directory. Be sure to add them manually to the job directory before running the job.\n\t{files_missing_str}"
            )

    def run(
        self,
        scheduler,
        job_dir_path=".",
        dock_files_dir_path=None,
        indock_file_path=None,
        custom_dock_executable=None,
        positives_tgz_file_path: Union[None, str] = None,
        negatives_tgz_file_path: Union[None, str] = None,
        retrodock_job_max_reattempts=0,
        retrodock_job_timeout_minutes: Union[None, str] = None,
        extra_submission_cmd_params_str: Union[None, str] = None,
        sleep_seconds_after_copying_output=0,
        export_negatives_mol2=False,
    ) -> None:
        """Run RetroDock job"""

        # validate args
        if dock_files_dir_path is None:
            dock_files_dir_path = os.path.join(job_dir_path, self.DOCK_FILES_DIR_NAME)
        if indock_file_path is None:
            indock_file_path = os.path.join(
                dock_files_dir_path, self.INDOCK_FILE_NAME
            )
        if positives_tgz_file_path is None:
            positives_tgz_file_path = os.path.join(
                job_dir_path, self.POSITIVES_TGZ_FILE_NAME
            )
        if negatives_tgz_file_path is None:
            negatives_tgz_file_path = os.path.join(job_dir_path, self.NEGATIVES_TGZ_FILE_NAME)
        if custom_dock_executable is None:
            dock_executable_path = DOCK3_EXECUTABLE_PATH
        else:
            dock_executable_path = custom_dock_executable

        try:
            File.validate_file_exists(positives_tgz_file_path)
            File.validate_file_exists(negatives_tgz_file_path)
        except FileNotFoundError:
            logger.error(
                "Positives TGZ file and/or negatives TGZ file not found. Did you put them in the job directory?\nNote: if you do not have positives and negatives, please use blastermaster instead of dockopt."
            )
            return
        if scheduler not in SCHEDULER_NAME_TO_CLASS_DICT:
            logger.error(
                f"scheduler flag must be one of: {list(SCHEDULER_NAME_TO_CLASS_DICT.keys())}"
            )
            return

        try:
            scheduler = SCHEDULER_NAME_TO_CLASS_DICT[scheduler]()
        except KeyError:
            logger.error(
                f"The following environmental variables are required to use the {scheduler} job scheduler: {SCHEDULER_NAME_TO_CLASS_DICT[scheduler].REQUIRED_ENV_VAR_NAMES}"
            )
            return

        try:
            TMPDIR = os.environ["TMPDIR"]
        except KeyError:
            logger.error(
                "The following environmental variables are required to submit retrodock jobs: TMPDIR"
            )
            return

        # create job directory
        job_dir = Dir(job_dir_path, create=True, reset=False)
        logger.info(f"Starting RetroDock job: {job_dir.path}")

        #
        retrospective_dataset = RetrospectiveDataset(positives_tgz_file_path, negatives_tgz_file_path, 'positives', 'negatives')

        dock_files_dir = Dir(dock_files_dir_path)
        dock_files = BlasterFiles(dock_files_dir).dock_files
        indock_file = IndockFile(indock_file_path)

        array_job_docking_configurations_file_path = os.path.join(job_dir.path, "array_job_docking_configurations.txt")
        with open(array_job_docking_configurations_file_path, 'w') as f:
            dockfile_paths_str = " ".join([dock_file.path for dock_file in astuple(dock_files)])
            f.write(f"{self.SINGLE_TASK_NUM} {indock_file.path} {dockfile_paths_str} {dock_executable_path}\n")

        #
        output_dir = Dir(os.path.join(job_dir.path, "output"), create=True, reset=False)

        #
        positives_dir = Dir(os.path.join(output_dir.path, "positives"), create=True, reset=False)
        positives_retrodock_job = ArrayDockingJob(
            name=f"retrodock_job_{uuid4()}_positives",
            job_dir=positives_dir,
            input_molecules_dir_path=retrospective_dataset.positives_dir_path,
            job_scheduler=scheduler,
            temp_storage_path=TMPDIR,
            array_job_docking_configurations_file_path=array_job_docking_configurations_file_path,
            job_timeout_minutes=retrodock_job_timeout_minutes,
            extra_submission_cmd_params_str=extra_submission_cmd_params_str,
            sleep_seconds_after_copying_output=sleep_seconds_after_copying_output,
            export_mol2=True,
        )

        #
        negatives_dir = Dir(os.path.join(output_dir.path, "negatives"), create=True, reset=False)
        negatives_retrodock_job = ArrayDockingJob(
            name=f"retrodock_job_{uuid4()}_negatives",
            job_dir=negatives_dir,
            input_molecules_dir_path=retrospective_dataset.negatives_dir_path,
            job_scheduler=scheduler,
            temp_storage_path=TMPDIR,
            array_job_docking_configurations_file_path=array_job_docking_configurations_file_path,
            job_timeout_minutes=retrodock_job_timeout_minutes,
            extra_submission_cmd_params_str=extra_submission_cmd_params_str,
            sleep_seconds_after_copying_output=sleep_seconds_after_copying_output,
            export_mol2=export_negatives_mol2,
        )

        #
        num_attempts_dict = collections.defaultdict(int)

        # submit jobs
        retrodock_jobs = [positives_retrodock_job, negatives_retrodock_job]
        for job in retrodock_jobs:
            sub_result, procs = job.submit_all_tasks(skip_if_complete=True)
            log_job_submission_result(job, sub_result, procs)
            if sub_result == JobSubmissionResult.SUCCESS:
                num_attempts_dict[job.name] += 1

        #
        def _resubmit_job_task(job: ArrayDockingJob):
            """Resubmit a job if it failed."""
            nonlocal num_attempts_dict
            if num_attempts_dict[job.name] < retrodock_job_max_reattempts + 1:
                sub_result, procs = job.submit_task(str(self.SINGLE_TASK_NUM), skip_if_complete=False)
                log_job_submission_result(job, sub_result, procs)
                num_attempts_dict[job.name] += 1
            else:
                raise Exception(f"Max job submission attempts ({retrodock_job_max_reattempts + 1}) exceeded. RetroDock job did not complete.")

        # wait for jobs to complete
        logger.info(f"Awaiting / processing tasks")
        while True:
            if any([job.is_on_job_scheduler_queue for job in retrodock_jobs]):
                time.sleep(5)
                continue
            else:
                if not all([job.is_complete for job in retrodock_jobs]):
                    for job in retrodock_jobs:
                        if job.task_failed(str(self.SINGLE_TASK_NUM)):
                            _resubmit_job_task(job)
                    continue

            #
            try:
                process_retrodock_job_results(
                    positives_outdock_file_path=os.path.join(positives_retrodock_job.job_dir.path, str(self.SINGLE_TASK_NUM), OUTDOCK_FILE_NAME),
                    negatives_outdock_file_path=os.path.join(negatives_retrodock_job.job_dir.path, str(self.SINGLE_TASK_NUM), OUTDOCK_FILE_NAME),
                    save_dir_path=job_dir.path,
                )
                logger.info(f"Successfully loaded both OUTDOCK files and processed results.")
            except Exception as e:
                logger.warning(f"Failed to parse output due to error: {e}")
                for job in retrodock_jobs:
                    _resubmit_job_task(job)

            #
            if all([job.is_complete for job in retrodock_jobs]):
                break

        #
        logger.info(f"Finished RetroDock job.")
