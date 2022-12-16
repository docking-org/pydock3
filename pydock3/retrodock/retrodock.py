import os
import logging
from uuid import uuid4
import time
from copy import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from joypy import joyplot

from pydock3.util import Script, CleanExit, get_dataclass_as_dict
from pydock3.files import (
    Dir,
    File,
    IndockFile,
    OutdockFile,
    INDOCK_FILE_NAME,
)
from pydock3.dockopt.roc import ROC
from pydock3.jobs import RetrodockJob, DOCK3_EXECUTABLE_PATH
from pydock3.blastermaster.blastermaster import BlasterFiles
from pydock3.jobs import JobSubmissionResult
from pydock3.job_schedulers import SGEJobScheduler, SlurmJobScheduler
from pydock3.retrodock import __file__ as RETRODOCK_INIT_FILE_PATH


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
SCHEDULER_NAME_TO_CLASS_DICT = {
    "sge": SGEJobScheduler,
    "slurm": SlurmJobScheduler,
}


# TODO: add this as a decorator of `submit` method of `DockoptJob`
def log_job_submission_result(job, submission_result, proc):
    if submission_result is JobSubmissionResult.SUCCESS:
        logger.info(f"Retrodock job '{job.name}' successfully submitted.\n")
    elif submission_result is JobSubmissionResult.FAILED:
        logger.info(f"Retrodock job submission failed for '{job.name}' due to error: {proc.stderr}\n")
    if submission_result is JobSubmissionResult.SKIPPED_BECAUSE_ALREADY_COMPLETE:
        logger.info(
            f"Retrodock job submission skipped for '{job.name}' since all its OUTDOCK files already exist.\n"
        )
    elif submission_result is JobSubmissionResult.SKIPPED_BECAUSE_STILL_RUNNING:
        logger.info(
            f"Retrodock job submission skipped for '{job.name}' since it is still running from a previous submission.\n"
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


def get_results_dataframe_from_actives_job_and_decoys_job_outdock_files(
    actives_outdock_file_path, decoys_outdock_file_path
):
    #
    actives_outdock_file = OutdockFile(actives_outdock_file_path)
    decoys_outdock_file = OutdockFile(decoys_outdock_file_path)

    #
    actives_outdock_df = actives_outdock_file.get_dataframe()
    decoys_outdock_df = decoys_outdock_file.get_dataframe()

    # set is_active column based on outdock file
    actives_outdock_df["is_active"] = [1 for _ in range(len(actives_outdock_df))]
    decoys_outdock_df["is_active"] = [0 for _ in range(len(decoys_outdock_df))]

    # set activity_class column based on outdock file
    actives_outdock_df["activity_class"] = [
        "active" for _ in range(len(actives_outdock_df))
    ]
    decoys_outdock_df["activity_class"] = [
        "decoy" for _ in range(len(decoys_outdock_df))
    ]

    # build dataframe of docking results from outdock files
    df = pd.DataFrame()
    df = pd.concat([df, actives_outdock_df], ignore_index=True)
    df = pd.concat([df, decoys_outdock_df], ignore_index=True)

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

    return df


class Retrodock(Script):
    JOB_DIR_NAME = "retrodock_job"
    DOCK_FILES_DIR_NAME = "dockfiles"
    INDOCK_FILE_NAME = "INDOCK"
    ACTIVES_TGZ_FILE_NAME = "actives.tgz"
    DECOYS_TGZ_FILE_NAME = "decoys.tgz"

    def __init__(self):
        super().__init__()

    def new(self, job_dir_path=JOB_DIR_NAME, overwrite=False):
        # create job dir
        job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create actives and decoys dirs
        actives_dir = Dir(path=job_dir_path, create=True, reset=False)
        decoys_dir = Dir(path=job_dir_path, create=True, reset=False)

    def run(
            self,
            scheduler,
            dock_executable_path=DOCK3_EXECUTABLE_PATH,
            job_dir_path=".",
            dock_files_dir_path=None,
            indock_file_path=None,
            actives_tgz_file_path=None,
            decoys_tgz_file_path=None,
            retrodock_job_max_reattempts=0,
            retrodock_job_timeout_minutes=None,
            max_scheduler_jobs_running_at_a_time=None,  # TODO
            export_decoy_poses=False,  # TODO
    ):
        # validate args
        if dock_files_dir_path is None:
            dock_files_dir_path = os.path.join(
                job_dir_path, self.DOCK_FILES_DIR_NAME
            )
        if indock_file_path is None:
            indock_file_path = os.path.join(
                dock_files_dir_path, self.INDOCK_FILE_NAME
            )  # TODO: come up with good placement of INDOCK file relative to dockfiles
        if actives_tgz_file_path is None:
            actives_tgz_file_path = os.path.join(
                job_dir_path, self.ACTIVES_TGZ_FILE_NAME
            )
        if decoys_tgz_file_path is None:
            decoys_tgz_file_path = os.path.join(job_dir_path, self.DECOYS_TGZ_FILE_NAME)
        try:
            File.validate_file_exists(actives_tgz_file_path)
            File.validate_file_exists(decoys_tgz_file_path)
        except FileNotFoundError:
            logger.error(
                "Actives TGZ file and/or decoys TGZ file not found. Did you put them in the job directory?\nNote: if you do not have actives and decoys, please use blastermaster instead of dockopt."
            )
            return
        if scheduler not in SCHEDULER_NAME_TO_CLASS_DICT:
            logger.error(
                f"scheduler flag must be one of: {list(SCHEDULER_NAME_TO_CLASS_DICT.keys())}"
            )
            return

        #
        try:
            scheduler = SCHEDULER_NAME_TO_CLASS_DICT[scheduler]()
        except KeyError:
            logger.error(
                f"The following environmental variables are required to use the {scheduler} job scheduler: {SCHEDULER_NAME_TO_CLASS_DICT[scheduler].REQUIRED_ENV_VAR_NAMES}"
            )
            return

        #
        try:
            TMPDIR = os.environ["TMPDIR"]
        except KeyError:
            logger.error(
                "The following environmental variables are required to submit retrodock jobs: TMPDIR"
            )
            return

        #
        job_dir = Dir(job_dir_path, create=True, reset=False)

        #
        if actives_tgz_file_path is not None:
            actives_tgz_file = File(path=actives_tgz_file_path)
        else:
            actives_tgz_file = None
        if actives_tgz_file_path is not None:
            decoys_tgz_file = File(path=decoys_tgz_file_path)
        else:
            decoys_tgz_file = None

        # write actives tgz and decoys tgz file paths to actives_and_decoys.sdi
        logger.info("Writing actives_and_decoys.sdi file...")
        retrodock_input_sdi_file = File(
            path=os.path.join(job_dir.path, "actives_and_decoys.sdi")
        )
        with open(retrodock_input_sdi_file.path, "w") as f:
            f.write(f"{actives_tgz_file.path}\n")
            f.write(f"{decoys_tgz_file.path}\n")
        logger.info("done")

        #
        dock_files = BlasterFiles(dock_files_dir_path).dock_files
        indock_file = IndockFile(indock_file_path)

        #
        retrodock_job_output_dir = Dir(
            path=os.path.join(job_dir.path, f"output"), create=True
        )

        #
        retrodock_job = RetrodockJob(
            name=f"retrodock_job_{uuid4()}",  # TODO: replace this with a hash based on docking configuration
            input_sdi_file=retrodock_input_sdi_file,
            dock_files=dock_files,
            indock_file=indock_file,
            output_dir=retrodock_job_output_dir,
            job_scheduler=scheduler,
            dock_executable_path=dock_executable_path,
            temp_storage_path=TMPDIR,
            max_reattempts=retrodock_job_max_reattempts,
        )
        sub_result, proc = retrodock_job.submit(
            job_timeout_minutes=retrodock_job_timeout_minutes,
            skip_if_complete=True,
        )
        log_job_submission_result(retrodock_job, sub_result, proc)

        #
        while not retrodock_job.is_complete:
            #
            if retrodock_job.is_running:
                time.sleep(1)
                continue

        #
        actives_outdock_file_path = os.path.join(
            retrodock_job.output_dir.path, "1", "OUTDOCK.0"
        )
        decoys_outdock_file_path = os.path.join(
            retrodock_job.output_dir.path, "2", "OUTDOCK.0"
        )

        # get dataframe of actives job results and decoys job results combined
        df = (
            get_results_dataframe_from_actives_job_and_decoys_job_outdock_files(
                actives_outdock_file_path, decoys_outdock_file_path
            )
        )

        #
        logger.info(
            f"Docking job '{retrodock_job.name}' completed. Successfully loaded OUTDOCK file(s)."
        )

        # sort dataframe by total energy score
        df = df.sort_values(
            by=["total_energy", "is_active"], na_position="last", ignore_index=True
        )  # sorting secondarily by 'is_active' (0 or 1) ensures that decoys are ranked before actives in case they have the same exact score (pessimistic approach)
        df = df.drop_duplicates(
            subset=["db2_file_path"], keep="first", ignore_index=True
        )  # keep only the best score per molecule

        # get ROC and calculate enrichment score of this job's docking set-up
        logger.debug("Calculating ROC and enrichment score...")
        booleans = df["is_active"]
        indices = df["total_energy"].fillna(
            np.inf
        )  # unscored molecules are assumed to have worst possible score (pessimistic approach)
        roc = ROC(booleans, indices)
        with open("enrichment_score", 'w') as f:
            f.write(f"{roc.enrichment_score}")
        logger.debug("done.")

        # write ROC plot image
        roc_plot_image_path = os.path.join(job_dir.path, "roc.png")
        roc.plot(save_path=roc_plot_image_path)

        # ridge plot for energy terms
        pivot_rows = []
        for i, row in df.iterrows():
            for col in [
                "total_energy",
                "electrostatic_energy",
                "vdw_energy",
                "polar_desolvation_energy",
                "apolar_desolvation_energy",
            ]:
                pivot_row = {"energy_term": col}
                if row["is_active"] == 1:
                    pivot_row["active"] = str_to_float(row[col])
                    pivot_row["decoy"] = np.nan
                else:
                    pivot_row["active"] = np.nan
                    pivot_row["decoy"] = str_to_float(row[col])
                pivot_rows.append(pivot_row)
        df_best_job_pivot = pd.DataFrame(pivot_rows)
        fig, ax = joyplot(
            data=df_best_job_pivot,
            by="energy_term",
            column=["active", "decoy"],
            color=["#686de0", "#eb4d4b"],
            legend=True,
            alpha=0.85,
            figsize=(12, 8),
            ylim="own",
        )
        plt.title("ridgeline plot: energy terms (actives vs. decoys)")
        plt.tight_layout()
        plt.savefig(os.path.join(job_dir.path, "energy.png"))
        plt.close(fig)

        # split violin plot of charges
        fig = plt.figure()
        sns.violinplot(
            data=df,
            x="charge",
            y="total_energy",
            split=True,
            hue="activity_class",
        )
        plt.title("split violin plot: charge (actives vs. decoys)")
        plt.tight_layout()
        plt.savefig(os.path.join(job_dir.path, "charge.png"))
        plt.close(fig)