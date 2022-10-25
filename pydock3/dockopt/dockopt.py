import itertools
import os
import shutil
import sys
from functools import wraps
from dataclasses import fields
from copy import copy
import logging
import collections
import time
import random
from datetime import datetime

import networkx as nx
import numpy as np
import pandas as pd
from dirhash import dirhash
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


from pydock3.util import (
    Script,
    CleanExit,
    get_dataclass_as_dict,
    validate_variable_type,
)
from pydock3.config import Parameter
from pydock3.blastermaster.blastermaster import BlasterFiles, get_blaster_steps
from pydock3.dockopt.config import DockoptParametersConfiguration
from pydock3.files import (
    Dir,
    File,
    IndockFile,
    OutdockFile,
    INDOCK_FILE_NAME,
)
from pydock3.blastermaster.util import (
    WorkingDir,
    BlasterFile,
    DockFiles,
    BlasterFileNames,
)
from pydock3.dockopt.roc import ROC
from pydock3.jobs import RetrodockJob, DOCK3_EXECUTABLE_PATH
from pydock3.job_schedulers import SlurmJobScheduler, SGEJobScheduler
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH
from pydock3.blastermaster.defaults import __file__ as DEFAULTS_INIT_FILE_PATH
from pydock3.blastermaster.programs.thinspheres.sph_lib import read_sph, write_sph


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#
plt.rcParams.update({"font.size": 14})

#
SCHEDULER_NAME_TO_CLASS_DICT = {
    "sge": SGEJobScheduler,
    "slurm": SlurmJobScheduler,
}

#
METRICS = ["enrichment_score"]
POSSIBLE_NON_PARAMETER_COLUMNS = METRICS + ["retrodock_job_num"]

#
RETRODOCK_JOB_DIR_COLUMN_NAME = "retrodock_job_dir"

#
ROC_IMAGE_FILE_NAME = "roc.png"


class Dockopt(Script):
    JOB_DIR_NAME = "dockopt_job"
    CONFIG_FILE_NAME = "dockopt_config.yaml"
    ACTIVES_TGZ_FILE_NAME = "actives.tgz"
    DECOYS_TGZ_FILE_NAME = "decoys.tgz"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(
        os.path.dirname(DOCKOPT_INIT_FILE_PATH), "default_dockopt_config.yaml"
    )
    DEFAULT_FILES_DIR_PATH = os.path.dirname(DEFAULTS_INIT_FILE_PATH)

    def __init__(self):
        super().__init__()

        self.job_dir = None  # assigned in .init()

    @staticmethod
    def handle_run_func(run_func):
        @wraps(run_func)
        def wrapper(self, *args, **kwargs):
            with CleanExit():
                logger.info(f"Running {self.__class__.__name__}")
                run_func(self, *args, **kwargs)

        return wrapper

    def init(self, job_dir_path=JOB_DIR_NAME, overwrite=False):
        # create job dir
        self.job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create working dir & copy in blaster files
        blaster_file_names = list(get_dataclass_as_dict(BlasterFileNames()).values())
        backup_blaster_file_paths = [
            os.path.join(self.DEFAULT_FILES_DIR_PATH, blaster_file_name)
            for blaster_file_name in blaster_file_names
        ]
        blaster_file_names_in_cwd = [f for f in blaster_file_names if os.path.isfile(f)]
        files_to_copy_str = "\n\t".join(blaster_file_names_in_cwd)
        if blaster_file_names_in_cwd:
            logger.info(
                f"Copying the following files from current directory into job working directory:\n\t{files_to_copy_str}"
            )
        else:
            logger.info(
                f"No blaster files detected in current working directory. Be sure to add them manually before running the job."
            )

        # copy in actives and decoys TGZ files
        tgz_files = [self.ACTIVES_TGZ_FILE_NAME, self.DECOYS_TGZ_FILE_NAME]
        tgz_file_names_in_cwd = [f for f in tgz_files if os.path.isfile(f)]
        tgz_file_names_not_in_cwd = [f for f in tgz_files if not os.path.isfile(f)]
        if tgz_file_names_in_cwd:
            files_to_copy_str = "\n\t".join(tgz_file_names_in_cwd)
            logger.info(
                f"Copying the following files from current directory into job directory:\n\t{files_to_copy_str}"
            )
            for tgz_file_name in tgz_file_names_in_cwd:
                self.job_dir.copy_in_file(tgz_file_name)
        if tgz_file_names_not_in_cwd:
            files_missing_str = "\n\t".join(tgz_file_names_not_in_cwd)
            logger.info(
                f"The following required files were not found in current working directory. Be sure to add them manually to the job directory before running the job.\n\t{files_missing_str}"
            )

        # write fresh config file from default file
        save_path = os.path.join(self.job_dir.path, self.CONFIG_FILE_NAME)
        DockoptParametersConfiguration.write_config_file(
            save_path, self.DEFAULT_CONFIG_FILE_PATH, overwrite=overwrite
        )

    @handle_run_func.__get__(0)
    def run(
        self,
        scheduler,
        job_dir_path=".",
        config_file_path=None,
        actives_tgz_file_path=None,
        decoys_tgz_file_path=None,
        retrodock_job_max_reattempts=0,
        retrodock_job_timeout_minutes=None,
        max_scheduler_jobs_running_at_a_time=None,  # TODO
        export_decoy_poses=False,  # TODO
    ):
        # validate args
        if config_file_path is None:
            config_file_path = os.path.join(job_dir_path, self.CONFIG_FILE_NAME)
        if actives_tgz_file_path is None:
            actives_tgz_file_path = os.path.join(
                job_dir_path, self.ACTIVES_TGZ_FILE_NAME
            )
        if decoys_tgz_file_path is None:
            decoys_tgz_file_path = os.path.join(job_dir_path, self.DECOYS_TGZ_FILE_NAME)
        try:
            File.validate_file_exists(config_file_path)
        except FileNotFoundError:
            logger.error("Config file not found. Are you in the job directory?")
            return
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
            temp_dir_path = os.environ["TMPDIR"]
        except KeyError:
            logger.error(
                "The following environmental variables are required to submit retrodock jobs: TMPDIR"
            )
            return

        #
        logger.info("Loading config file...")
        config = DockoptParametersConfiguration(config_file_path)  # TODO
        logger.info("done.")

        #
        config_params_str = "\n".join(  # TODO
            [
                f"{param_name}: {param.value}"
                for param_name, param in config.param_dict.items()
            ]
        )
        logger.info(f"Parameters:\n{config_params_str}")

        # TODO


class Criterion(object):
    def __init__(self):
        pass

    @property
    def name(self):
        raise NotImplementedError

    def calculate(self, runnable):
        raise NotImplementedError


class EnrichmentScore(Criterion):
    def __init__(self):
        super().__init__()

    @property
    def name(self):
        return "enrichment_score"

    def calculate(self, runnable):
        # TODO
        pass


def record_started_utc_and_finished_utc_for_run_method(_cls):
    run = getattr(_cls, "run")

    def new_run(self):
        self.started_utc = datetime.utcnow()
        run(self)
        self.finished_utc = datetime.utcnow()

    setattr(_cls, "run", new_run)

    return _cls


class RunnableJob(object):
    def __init__(self, job_dir_path, criterion, top_n_jobs_to_keep):
        #
        self.criterion = criterion
        self.top_n_job_to_keep = top_n_jobs_to_keep

        #
        self.job_dir = Dir(job_dir_path)
        self.results_csv_file_path = os.path.join(self.job_dir.path, "results.csv")

        #
        self.started_utc = None  # set by record_started_utc_and_finished_utc_for_run_method() decorator; see .__init_subclass__()
        self.finished_utc = None  # set by record_started_utc_and_finished_utc_for_run_method() decorator; see .__init_subclass__()

    def __init_subclass__(cls, **kwargs):
        return record_started_utc_and_finished_utc_for_run_method(_cls=cls)

    def run(self, *args, **kwargs):
        raise NotImplementedError

    def write_results_csv(self, *args, **kwargs):
        raise NotImplementedError

    def load_results_csv(self):
        df = pd.read_csv(self.results_csv_file_path)
        df = df.loc[
            :, ~df.columns.str.contains("^Unnamed")
        ]  # remove useless index column
        df = df.sort_values(by=self.criterion.name, ascending=False, ignore_index=True)

        return df

    def get_n_best_jobs_dir_paths(self):
        df = self.load_results_csv()

        return df.nlargest(self.top_n_job_to_keep, self.criterion.name)[
            RETRODOCK_JOB_DIR_COLUMN_NAME
        ].tolist()


class RunnableJobWithReport(RunnableJob):
    def __init__(self, job_dir_path, criterion, top_n_jobs_to_keep):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
        )

        #
        self.report_pdf_file_path = os.path.join(self.job_dir.path, "report.pdf")

    def run(self, *args, **kwargs):
        raise NotImplementedError

    def write_results_csv(self, *args, **kwargs):
        raise NotImplementedError

    def write_report_pdf(self):
        def get_ordinal(n):
            return "%d%s" % (
                n,
                "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10 :: 4],
            )

        df = self.load_results_csv()

        with PdfPages(self.report_pdf_file_path) as f:
            #
            fig = plt.figure(figsize=(11.0, 8.5))

            # write box plot for data grouped by runnable sequence identifier
            # TODO

            #
            for i, best_job_dir_path in enumerate(self.get_n_best_jobs_dir_paths()):
                if File.file_exists(
                    os.path.join(best_job_dir_path, ROC_IMAGE_FILE_NAME)
                ):
                    image = mpimg.imread(f)
                    plt.axis("off")
                    plt.suptitle(
                        f"linear-log ROC plot of {get_ordinal(i+1)} best job\n{best_job_dir_path}"
                    )
                    plt.imshow(image)
            f.savefig(fig, bbox_inches="tight")
            plt.close(fig)

            #
            multivalued_config_param_columns = [
                column
                for column in df.columns
                if column not in POSSIBLE_NON_PARAMETER_COLUMNS
                and df[column].nunique() > 1
            ]
            for column in multivalued_config_param_columns:
                fig = plt.figure(figsize=(8, 6))

                if len(df) == len(df.drop_duplicates(subset=[column])):
                    df.plot.bar(x=column, y=self.criterion.name, rot=0)
                else:
                    sns.boxplot(
                        data=df,
                        x=column,
                        y=self.criterion.name,
                        showfliers=False,
                        boxprops={"facecolor": "None"},
                    )
                    sns.stripplot(data=df, x=column, y=self.criterion.name, zorder=0.5)
                    fig.autofmt_xdate(rotation=25)
                    plt.yticks(rotation=0)

                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)

            #
            df = df.sort_values(
                by=self.criterion.name, ascending=False, ignore_index=True
            )
            for column_1, column_2 in itertools.combinations(
                multivalued_config_param_columns, 2
            ):
                #
                fig, ax = plt.subplots()
                fig.set_size_inches(8.0, 6.0)

                #
                df_no_duplicates = df.drop_duplicates(
                    subset=[column_1, column_2], keep="first", ignore_index=True
                )

                #
                df_pivot = pd.pivot_table(
                    df_no_duplicates,
                    values=self.criterion.name,
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
                    cmap="icefire",
                    robust=True,
                    cbar_kws={"label": self.criterion.name},
                )
                fig.autofmt_xdate(rotation=25)
                plt.yticks(rotation=0)

                f.savefig(fig, bbox_inches="tight")
                plt.close(fig)


class DockingConfigurationNodeBranchingJob(RunnableJobWithReport):
    WORKING_DIR_NAME = "working"
    RETRODOCK_JOBS_DIR_NAME = "retrodock_jobs"
    BEST_RETRODOCK_JOBS_DIR_NAME = "best_retrodock_jobs"

    def __init__(
        self, job_dir_path, criterion, top_n_jobs_to_keep, blaster_files_to_copy_in
    ):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
        )

        #
        backup_blaster_files = []  # TODO
        self.working_dir = WorkingDir(
            path=os.path.join(self.job_dir.path, self.WORKING_DIR_NAME),
            create=True,
            reset=False,
            files_to_copy_in=blaster_files_to_copy_in,
            backup_files_to_copy_in=backup_blaster_files,
        )
        self.retrodock_jobs_dir = Dir(
            path=os.path.join(self.job_dir.path, self.RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=False,
        )
        self.best_retrodock_jobs_dir = Dir(
            path=os.path.join(self.job_dir.path, self.BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.config_yaml_dict = None  # TODO

        #
        self.actives_tgz_file = None  # set at beginning of .run()
        self.decoys_tgz_file = None  # set at beginning of .run()

        #
        dock_files_generation_param_dicts = None  # TODO

        #
        self.blaster_files = BlasterFiles(working_dir=self.working_dir)

        # get directed acyclical graph defining how to get all combinations of dock files we need from the provided input files and parameters
        graph = nx.DiGraph()
        dock_file_nodes_combinations = []
        for dock_files_generation_param_dict in dock_files_generation_param_dicts:
            # get param_dict for get_blaster_steps
            # each value in dict must be an instance of Parameter
            def flatten_dict_into_parameter_dict(d, key_prefix=""):
                new_d = {}
                for key, value in d.items():
                    this_key = f"{key_prefix}{key}"
                    if isinstance(value, dict):
                        new_d.update(
                            flatten_dict_into_parameter_dict(value, f"{this_key}.")
                        )
                    else:
                        new_d[this_key] = Parameter(this_key, value)
                return new_d

            #
            steps = get_blaster_steps(
                blaster_files=self.blaster_files,
                param_dict=flatten_dict_into_parameter_dict(
                    dock_files_generation_param_dict
                ),
                working_dir=self.working_dir,
            )

            # form subgraph for this dock_files_generation_param_dict from the blaster steps it defines
            subgraph = nx.DiGraph()
            for step in steps:

                # add infile nodes
                for infile in step.infiles:
                    if (
                        self._get_blaster_file_node_with_same_name(
                            infile.original_file_in_working_dir.name, subgraph
                        )
                        is not None
                    ):
                        continue
                    subgraph.add_node(
                        self._get_hash_for_new_blaster_file_node(
                            infile, step, subgraph
                        ),
                        blaster_file=infile.original_file_in_working_dir,
                    )

                # add outfile nodes
                for outfile in step.outfiles:
                    if self._get_blaster_file_node_with_same_name(
                        outfile.original_file_in_working_dir.name, subgraph
                    ):
                        raise Exception(
                            f"Attempting to add outfile to subgraph that already has said outfile as node: {outfile.original_file_in_working_dir.name}"
                        )
                    subgraph.add_node(
                        self._get_hash_for_new_blaster_file_node(
                            outfile, step, subgraph
                        ),
                        blaster_file=outfile.original_file_in_working_dir,
                    )

                # add parameter nodes
                for parameter in step.parameters:
                    subgraph.add_node(parameter.__hash__(), parameter=parameter)

                # get step hash
                infiles_dict_items_list = step.infiles._asdict().items()
                outfiles_dict_items_list = step.outfiles._asdict().items()
                parameters_dict_items_list = step.parameters._asdict().items()
                step_hash = hash(
                    tuple(
                        sorted(
                            [
                                (
                                    infile_step_var_name,
                                    self._get_blaster_file_node_with_same_name(
                                        infile.original_file_in_working_dir.name,
                                        subgraph,
                                    ),
                                )
                                for infile_step_var_name, infile in infiles_dict_items_list
                            ]
                            + [
                                (
                                    outfile_step_var_name,
                                    self._get_blaster_file_node_with_same_name(
                                        outfile.original_file_in_working_dir.name,
                                        subgraph,
                                    ),
                                )
                                for outfile_step_var_name, outfile in outfiles_dict_items_list
                            ]
                            + [
                                (parameter_step_var_name, parameter.__hash__())
                                for parameter_step_var_name, parameter in parameters_dict_items_list
                            ]
                        )
                    )
                )

                # connect each infile node to every outfile node
                for (infile_step_var_name, infile), (
                    outfile_step_var_name,
                    outfile,
                ) in itertools.product(
                    infiles_dict_items_list, outfiles_dict_items_list
                ):
                    subgraph.add_edge(
                        infile.original_file_in_working_dir.name,
                        outfile.original_file_in_working_dir.name,
                        step_class=step.__class__,
                        step_hash=step_hash,
                        parent_node_step_var_name=infile_step_var_name,
                        child_node_step_var_name=outfile_step_var_name,
                    )

                # connect each parameter node to every outfile nodes
                for (parameter_step_var_name, parameter), (
                    outfile_step_var_name,
                    outfile,
                ) in itertools.product(
                    parameters_dict_items_list, outfiles_dict_items_list
                ):
                    subgraph.add_edge(
                        parameter.name,
                        outfile.original_file_in_working_dir.name,
                        step_class=step.__class__,
                        step_hash=step_hash,
                        parent_node_step_var_name=parameter_step_var_name,
                        child_node_step_var_name=outfile_step_var_name,
                    )

                # record the combination of dock file nodes for this subgraph
                dock_file_nodes_combinations.append(
                    sorted(
                        (
                            [
                                node
                                for node in subgraph.nodes
                                if subgraph.out_degree(node) == 0
                            ]
                        )
                    )
                )

                # merge subgraph into full graph
                graph = nx.compose(graph, subgraph)

        #
        step_hash_to_edges_dict = collections.defaultdict(list)
        step_hash_to_step_class_dict = {}
        step_hash_to_step_class_instance_dict = {}
        for u, v, data in graph.edges(data=True):
            step_hash_to_edges_dict[data["step_hash"]].append((u, v))
            step_hash_to_step_class_dict[data["step_hash"]] = data["step_class"]

        #
        step_class_name_to_num_unique_instances_witnessed_so_far_counter = (
            collections.defaultdict(int)
        )
        step_hash_to_step_dir_path_dict = {}
        for step_hash, edges in step_hash_to_edges_dict.items():
            step_dir_path = graph.get_edge_data(*edges[0])["step_dir"].path
            step_class_name = graph.get_edge_data(*edges[0])["step_class"].__name__
            if (
                step_class_name
                not in step_class_name_to_num_unique_instances_witnessed_so_far_counter
            ):
                step_class_name_to_num_unique_instances_witnessed_so_far_counter[
                    step_class_name
                ] += 1
            if step_hash not in step_hash_to_step_dir_path_dict:
                step_hash_to_step_dir_path_dict[step_hash] = os.path.join(
                    f"{step_dir_path}_{step_class_name_to_num_unique_instances_witnessed_so_far_counter[step_class_name]}"
                )

        #
        for step_hash, edges in step_hash_to_edges_dict.items():
            #
            step_class = graph.get_edge_data(*edges[0])["step_class"]

            #
            step_dir = Dir(path=step_hash_to_step_dir_path_dict[step_hash])

            #
            kwargs = {"step_dir": step_dir}
            for (parent_node, child_node) in edges:
                graph[parent_node][child_node][
                    "step_instance"
                ] = step_hash_to_step_class_instance_dict[step_hash]

                edge_data_dict = graph.get_edge_data(parent_node, child_node)
                parent_node_data_dict = graph.nodes[parent_node]
                child_node_data_dict = graph.nodes[child_node]
                parent_node_step_var_name = edge_data_dict["parent_node_step_var_name"]
                child_node_step_var_name = edge_data_dict["child_node_step_var_name"]
                if "blaster_file" in parent_node_data_dict:
                    kwargs[parent_node_step_var_name] = parent_node_data_dict[
                        "blaster_file"
                    ]
                if "parameter" in parent_node_data_dict:
                    kwargs[parent_node_step_var_name] = parent_node_data_dict[
                        "parameter"
                    ]
                if "blaster_file" in child_node_data_dict:
                    kwargs[child_node_step_var_name] = child_node_data_dict[
                        "blaster_file"
                    ]
                if "parameter" in child_node_data_dict:
                    kwargs[child_node_step_var_name] = child_node_data_dict["parameter"]

            #
            step_hash_to_step_class_instance_dict[step_hash] = step_class(**kwargs)

            #
            for (parent_node, child_node) in edges:
                graph.get_edge_data(parent_node, child_node)[
                    "step_instance"
                ] = step_hash_to_step_class_instance_dict[step_hash]

        # validate that there are no cycles (i.e. that it is a DAG)
        if not nx.is_directed_acyclic_graph(graph):
            raise Exception("Cycle found in blaster targets DAG!")

        #
        self.graph = graph
        logger.debug(
            f"Graph initialized with:\n\tNodes: {self.graph.nodes}\n\tEdges: {self.graph.edges}"
        )

        #
        self.dock_file_nodes_combinations = dock_file_nodes_combinations

    def run(
        self,
        scheduler,
        temp_dir_path,
        actives_tgz_file_path=None,
        decoys_tgz_file_path=None,
        retrodock_job_max_reattempts=0,
        retrodock_job_timeout_minutes=None,
        max_scheduler_jobs_running_at_a_time=None,  # TODO
        export_decoy_poses=False,  # TODO
    ):
        #
        if actives_tgz_file_path is not None:
            self.actives_tgz_file = File(path=actives_tgz_file_path)
        else:
            self.actives_tgz_file = None
        if actives_tgz_file_path is not None:
            self.decoys_tgz_file = File(path=decoys_tgz_file_path)
        else:
            self.decoys_tgz_file = None

        # run edge steps
        logger.info("Running blaster steps / validating blaster files")
        for edge in nx.edge_bfs(self.graph, self._get_start_nodes(self.graph)):
            self._run_edge_step(edge, self.graph)

        # matching spheres perturbation
        if self.param_dict["matching_spheres_perturbation.use"].value:
            #
            unperturbed_file_name_to_perturbed_file_names_dict = (
                collections.defaultdict(list)
            )
            blaster_file_nodes = [
                node_name
                for node_name, node_data in self.graph.nodes(data=True)
                if node_data.get("blaster_file")
            ]
            for node_name, node_data in self.graph.nodes(data=True):
                if (
                    node_data["blaster_file"].name
                    == self.blaster_files.matching_spheres_file.name
                ):
                    spheres = read_sph(
                        os.path.join(self.working_dir.path, node_name),
                        chosen_cluster="A",
                        color="A",
                    )

                    #
                    for i in range(
                        int(
                            self.param_dict[
                                "matching_spheres_perturbation.num_samples_per_matching_spheres_file"
                            ].value
                        )
                    ):
                        #
                        perturbed_file_name = f"{node_name}_{i + 1}"
                        perturbed_file_path = os.path.join(
                            self.working_dir.path, perturbed_file_name
                        )
                        unperturbed_file_name_to_perturbed_file_names_dict[
                            node_name
                        ].append(perturbed_file_name)

                        # skip perturbation if perturbed file already exists
                        if File.file_exists(perturbed_file_path):
                            continue

                        # perturb all spheres in file
                        new_spheres = []
                        for sphere in spheres:
                            new_sphere = copy(sphere)
                            max_deviation = float(
                                config.param_dict[
                                    "matching_spheres_perturbation.max_deviation_angstroms"
                                ].value
                            )
                            perturbation_xyz = tuple(
                                [
                                    random.uniform(
                                        -max_deviation,
                                        max_deviation,
                                    )
                                    for _ in range(3)
                                ]
                            )
                            new_sphere.X += perturbation_xyz[0]
                            new_sphere.Y += perturbation_xyz[1]
                            new_sphere.Z += perturbation_xyz[2]
                            new_spheres.append(new_sphere)

                        # write perturbed spheres to new matching spheres file
                        write_sph(perturbed_file_path, new_spheres)

            #
            dock_files_combinations_after_modifications = []
            dock_files_generation_param_dicts_after_modifications = []
            for i, (
                dock_files_combination,
                dock_files_generation_param_dict,
            ) in enumerate(
                zip(
                    self.dock_files_combinations, self.dock_files_generation_param_dicts
                )
            ):
                for (
                    perturbed_file_name
                ) in unperturbed_file_name_to_perturbed_file_names_dict[
                    dock_files_combination.matching_spheres_file.name
                ]:
                    new_dock_files_combination = copy(dock_files_combination)
                    new_dock_files_combination.matching_spheres_file = BlasterFile(
                        os.path.join(self.working_dir.path, perturbed_file_name)
                    )
                    dock_files_combinations_after_modifications.append(
                        new_dock_files_combination
                    )
                    dock_files_generation_param_dicts_after_modifications.append(
                        self.dock_files_generation_param_dicts[i]
                    )
        else:
            dock_files_combinations_after_modifications = self.dock_files_combinations
            dock_files_generation_param_dicts_after_modifications = (
                self.dock_files_generation_param_dicts
            )

        #
        matching_spheres_perturbation_param_dict = {  # TODO
            key: value
            for key, value in self.param_dict.items()
            if key.startswith("matching_spheres_perturbation.")
        }

        #
        # TODO: there is an assumption here that none of the indock config params are used in any of the blaster steps; if this assumption ever breaks, this method will need to be replaced.
        job_param_dicts_indock_subset = [
            {
                key: value
                for key, value in job_param_dict.items()
                if key.startswith("indock.")
            }
            for job_param_dict in config.job_param_dicts
        ]
        job_param_dicts_indock_subset = [
            dict(s)
            for s in set(
                frozenset(job_param_dict.items())
                for job_param_dict in job_param_dicts_indock_subset
            )
        ]  # get unique dicts

        #
        if isinstance(config.param_dict["custom_dock_executable"].value, list):
            dock_executable_paths = []
            for dock_executable_path in config.param_dict[
                "custom_dock_executable"
            ].value:
                if dock_executable_path is None:
                    dock_executable_paths.append(DOCK3_EXECUTABLE_PATH)
                else:
                    dock_executable_paths.append(dock_executable_path)
        else:
            if config.param_dict["custom_dock_executable"].value is None:
                dock_executable_paths = [DOCK3_EXECUTABLE_PATH]
            else:
                dock_executable_paths = [
                    config.param_dict["custom_dock_executable"].value
                ]

        #
        docking_configuration_info_combinations = list(
            itertools.product(
                dock_executable_paths,
                zip(
                    dock_files_combinations_after_modifications,
                    dock_files_generation_param_dicts_after_modifications,
                ),
                job_param_dicts_indock_subset,
            )
        )

        # make indock file for each combination of (1) set of dock files and (2) job_param_dict_indock_subset
        logger.info("Making INDOCK files...")
        parameter_dicts = []
        docking_configurations = []
        for i, (
            dock_executable_path,
            (dock_files, dock_files_generation_param_dict),
            job_param_dict_indock_subset,
        ) in enumerate(docking_configuration_info_combinations):
            # get full parameter dict
            parameter_dict = {p.name: p.value for p in dock_files_generation_param_dict}
            parameter_dict.update(
                matching_spheres_perturbation_param_dict
            )  # add matching spheres perturbation params
            parameter_dict.update(job_param_dict_indock_subset)  # add indock params
            parameter_dict["dock_executable_path"] = dock_executable_path

            # make indock file for each combination of dock files
            indock_file_name = f"{INDOCK_FILE_NAME}_{i + 1}"
            indock_file = IndockFile(
                path=os.path.join(self.working_dir.path, indock_file_name)
            )
            indock_file.write(dock_files, parameter_dict)

            #
            parameter_dicts.append(parameter_dict)
            docking_configurations.append(
                (dock_executable_path, dock_files, indock_file)
            )

        #
        all_docking_configuration_file_names = []
        for dock_executable_path, dock_files, indock_file in docking_configurations:
            all_docking_configuration_file_names.append(dock_executable_path)
            dock_file_names = [
                getattr(dock_files, dock_file_field.name).name
                for dock_file_field in fields(dock_files)
            ]
            all_docking_configuration_file_names += dock_file_names
            all_docking_configuration_file_names.append(indock_file.name)
        all_docking_configuration_file_names = list(
            set(all_docking_configuration_file_names)
        )

        #
        job_hash = dirhash(
            self.working_dir.path,
            "md5",
            match=all_docking_configuration_file_names,
        )
        with open(os.path.join(self.job_dir.path, "job_hash.md5"), "w") as f:
            f.write(f"{job_hash}\n")

        # write actives tgz and decoys tgz file paths to actives_and_decoys.sdi
        logger.info("Writing actives_and_decoys.sdi file...")
        retrodock_input_sdi_file = File(
            path=os.path.join(self.job_dir.path, "actives_and_decoys.sdi")
        )
        with open(retrodock_input_sdi_file.path, "w") as f:
            f.write(f"{self.actives_tgz_file.path}\n")
            f.write(f"{self.decoys_tgz_file.path}\n")
        logger.info("done")

        #
        retrodock_jobs = []
        retrodock_job_dirs = []
        retrodock_job_num_to_docking_configuration_file_names_dict = {}
        for i, (dock_executable_path, dock_files, indock_file) in enumerate(
            docking_configurations
        ):
            #
            retro_dock_job_num = str(i + 1)
            docking_configuration_file_names = [
                getattr(dock_files, dock_file_field.name).name
                for dock_file_field in fields(dock_files)
            ] + [indock_file.name]
            retrodock_job_num_to_docking_configuration_file_names_dict[
                retro_dock_job_num
            ] = docking_configuration_file_names

            #
            retrodock_job_dir = Dir(
                path=os.path.join(self.retrodock_jobs_dir.path, retro_dock_job_num),
                create=True,
            )
            retrodock_job_output_dir = Dir(
                path=os.path.join(retrodock_job_dir.path, f"output"), create=True
            )

            #
            retrodock_job = RetrodockJob(
                name=f"dockopt_job_{job_hash}_{retrodock_job_dir.name}",
                input_sdi_file=retrodock_input_sdi_file,
                dock_files=dock_files,
                indock_file=indock_file,
                output_dir=retrodock_job_output_dir,
                job_scheduler=scheduler,
                dock_executable_path=dock_executable_path,
                temp_storage_path=temp_dir_path,
                max_reattempts=retrodock_job_max_reattempts,
            )
            retrodock_jobs.append(retrodock_job)
            retrodock_job_dirs.append(retrodock_job_dir)
        logger.debug("done")

        #
        def submit_retrodock_job(retrodock_job, skip_if_complete):
            logger.info(f"Submitting docking job for {retrodock_job.name}...")
            proc = retrodock_job.run(
                job_timeout_minutes=retrodock_job_timeout_minutes,
                skip_if_complete=skip_if_complete,
            )
            if proc is None:
                logger.info(
                    f"Skipping docking job submission for {retrodock_job.name} since all its OUTDOCK files already exist.\n"
                )
            else:
                logger.debug(
                    f"Retrodock job submission system call returned: {proc}\n\nstdout:{proc.stdout}\n\nstderr:{proc.stderr}\n"
                )
                if proc.stderr:
                    logger.info(f"Job submission failed due to error: {proc.stderr}\n")
                else:
                    logger.info("done.\n")

        # submit docking jobs
        for retrodock_job in retrodock_jobs:
            submit_retrodock_job(retrodock_job, skip_if_complete=True)

        # make a queue of tuples containing job-relevant data for processing
        RetrodockJobInfoTuple = collections.namedtuple(
            "RetrodockJobInfoTuple", "job job_dir parameter_dict"
        )
        retrodock_jobs_processing_queue = [
            RetrodockJobInfoTuple(
                retrodock_jobs[i], retrodock_job_dirs[i], parameter_dicts[i]
            )
            for i in range(len(retrodock_jobs))
        ]

        # process results of docking jobs
        logger.info(
            f"Awaiting / processing retrodock job results ({len(retrodock_job_dirs)} jobs in total)"
        )
        data_dicts = []
        while len(retrodock_jobs_processing_queue) > 0:
            #
            retrodock_job_info_tuple = retrodock_jobs_processing_queue.pop(0)
            retrodock_job, retrodock_job_dir, parameter_dict = retrodock_job_info_tuple

            #
            if retrodock_job.is_running:
                retrodock_jobs_processing_queue.append(
                    retrodock_job_info_tuple
                )  # move job to back of queue
                time.sleep(
                    1
                )  # sleep a bit while waiting for outdock file in order to avoid wasteful queue-cycling
                continue  # move on to next job in queue while job continues to run
            else:
                if (
                    not retrodock_job.is_complete
                ):  # not all expected OUTDOCK files exist yet
                    time.sleep(
                        1
                    )  # sleep for a bit and check again in case job just finished
                    if not retrodock_job.is_complete:
                        # job must have timed out / failed
                        logger.warning(
                            f"Job failure / time out witnessed for job: {retrodock_job.name}"
                        )
                        if retrodock_job.num_attempts > retrodock_job_max_reattempts:
                            logger.warning(
                                f"Max job reattempts exhausted for job: {retrodock_job.name}"
                            )
                            continue  # move on to next job in queue without re-attempting failed job
                        submit_retrodock_job(
                            retrodock_job, skip_if_complete=False
                        )  # re-attempt job
                        retrodock_jobs_processing_queue.append(
                            retrodock_job_info_tuple
                        )  # move job to back of queue
                        continue  # move on to next job in queue while docking job runs

            #
            actives_outdock_file_path = os.path.join(
                retrodock_job.output_dir.path, "1", "OUTDOCK.0"
            )
            decoys_outdock_file_path = os.path.join(
                retrodock_job.output_dir.path, "2", "OUTDOCK.0"
            )

            # load outdock file and get dataframe
            actives_outdock_file = OutdockFile(actives_outdock_file_path)
            decoys_outdock_file = OutdockFile(decoys_outdock_file_path)
            try:
                actives_outdock_df = actives_outdock_file.get_dataframe()
                decoys_outdock_df = decoys_outdock_file.get_dataframe()
            except Exception as e:  # if outdock file failed to be parsed then re-attempt job
                logger.warning(f"Failed to parse outdock file(s) due to error: {e}")
                if retrodock_job.num_attempts > retrodock_job_max_reattempts:
                    logger.warning(
                        f"Max job reattempts exhausted for job: {retrodock_job.name}"
                    )
                    continue  # move on to next job in queue without re-attempting failed job
                submit_retrodock_job(
                    retrodock_job, skip_if_complete=False
                )  # re-attempt job
                retrodock_jobs_processing_queue.append(
                    retrodock_job_info_tuple
                )  # move job to back of queue
                continue  # move on to next job in queue while docking job runs

            #
            logger.info(
                f"Docking job '{retrodock_job.name}' completed. Successfully loaded OUTDOCK file(s)."
            )

            # set is_active column based on outdock file
            actives_outdock_df["is_active"] = [
                1 for _ in range(len(actives_outdock_df))
            ]
            decoys_outdock_df["is_active"] = [0 for _ in range(len(decoys_outdock_df))]

            # build dataframe of docking results from outdock files
            df = pd.DataFrame()
            df = pd.concat([df, actives_outdock_df], ignore_index=True)
            df = pd.concat([df, decoys_outdock_df], ignore_index=True)

            # sort dataframe by total energy score
            df["Total"] = df["Total"].astype(float)
            df = df.sort_values(
                by=["Total", "is_active"], na_position="last", ignore_index=True
            )  # sorting secondarily by 'is_active' (0 or 1) ensures that decoys are ranked before actives in case they have the same exact score (pessimistic approach)
            df = df.drop_duplicates(
                subset=["db2_file_path"], keep="first", ignore_index=True
            )

            # make data dict for this job (will be used to make dataframe for results of all jobs)
            data_dict = copy(parameter_dict)
            data_dict["retrodock_job_num"] = retrodock_job_dir.name

            # get ROC and calculate enrichment score of this job's docking set-up
            logger.debug("Calculating ROC and enrichment score...")
            booleans = df["is_active"]
            indices = df["Total"].fillna(
                np.inf
            )  # unscored molecules are assumed to have worst possible score (pessimistic approach)
            roc = ROC(booleans, indices)
            data_dict["enrichment_score"] = roc.enrichment_score
            logger.debug("done.")

            # write ROC plot image
            roc_plot_image_path = os.path.join(
                retrodock_job_dir.path, ROC_IMAGE_FILE_NAME
            )
            roc.plot(save_path=roc_plot_image_path)

            # save data_dict for this job
            data_dicts.append(data_dict)

        # write jobs completion status
        num_jobs_completed = len(
            [1 for retrodock_job in retrodock_jobs if retrodock_job.is_complete]
        )
        logger.info(
            f"Finished {num_jobs_completed} out of {len(retrodock_jobs)} retrodock jobs."
        )
        if num_jobs_completed == 0:
            logger.error(
                "All retrodock jobs failed. Something is wrong. Please check logs."
            )
            return

        # make dataframe of optimization job results
        df = pd.DataFrame(data=data_dicts)
        df = df.sort_values(by=self.criterion.name, ascending=False, ignore_index=True)

        # save optimization job results dataframe to csv
        optimization_results_csv_file_path = os.path.join(
            self.job_dir.path, "dockopt_job_results.csv"
        )
        logger.debug(
            f"Saving optimization job results to {optimization_results_csv_file_path}"
        )
        df.to_csv(optimization_results_csv_file_path)

        # copy best job to output dir
        best_retrodock_job_dir = Dir(
            os.path.join(self.job_dir.path, "best_retrodock_job")
        )  # TODO
        logger.debug(
            f"Copying dockfiles of best job results to {best_retrodock_job_dir.path}"
        )
        if os.path.isdir(best_retrodock_job_dir.path):
            shutil.rmtree(best_retrodock_job_dir.path, ignore_errors=True)
        best_retrodock_job_num = df["retrodock_job_num"]
        shutil.copytree(
            os.path.join(self.retrodock_jobs_dir.path, best_retrodock_job_num),
            best_retrodock_job_dir.path,
        )

        # copy docking configuration files to best jobs dir
        best_retrodock_job_dockfiles_dir = Dir(  # TODO
            os.path.join(best_retrodock_job_dir.path, "dockfiles"), create=True
        )
        for file_name in retrodock_job_num_to_docking_configuration_file_names_dict[
            best_retrodock_job_num
        ]:
            best_retrodock_job_dockfiles_dir.copy_in_file(
                os.path.join(self.working_dir.path, file_name)
            )

        # write report PDF
        self.write_report_pdf()

    def write_results_csv(self, df):
        df.to_csv()

    @staticmethod
    def _get_blaster_file_node_with_same_name(name, g):
        blaster_file_nodes = [
            g.nodes[node] for node in g.nodes if g.nodes[node].get("blaster_file")
        ]
        if len(blaster_file_nodes) == 0:
            return None
        blaster_file_nodes_with_same_name = [
            blaster_file_node
            for blaster_file_node in blaster_file_nodes
            if name == blaster_file_node["blaster_file"].name
        ]
        if len(blaster_file_nodes_with_same_name) == 0:
            return None
        (blaster_file_node_with_same_name,) = blaster_file_nodes_with_same_name

        return blaster_file_node_with_same_name

    @staticmethod
    def _get_hash_for_new_blaster_file_node(blaster_file, step, g):
        blaster_file_node_with_same_name = (
            DockingConfigurationNodeBranchingJob._get_blaster_file_node_with_same_name(
                blaster_file.original_file_in_working_dir.name, g
            )
        )
        if blaster_file_node_with_same_name is None:
            return hash(blaster_file.original_file_in_working_dir.name)
        else:
            parent_nodes = [
                DockingConfigurationNodeBranchingJob._get_blaster_file_node_with_same_name(
                    infile.original_file_in_working_dir.name, g
                )
                for infile in step.infiles
            ]
            return hash(tuple(sorted([step.__class__.__name__] + parent_nodes)))

    @staticmethod
    def _get_start_nodes(g):
        start_nodes = []
        for node in g.nodes:
            if g.in_degree(node) == 0:
                start_nodes.append(node)
        return start_nodes

    @staticmethod
    def _run_edge_step(edge, g):
        if g.get_edge_data(*edge).get("step").is_done:
            return

        _, child_node = edge
        all_parent_nodes = g.predecessors(child_node)
        for parent_node in all_parent_nodes:
            parent_parent_nodes = g.predecessors(parent_node)
            for parent_parent_node in parent_parent_nodes:
                DockingConfigurationNodeBranchingJob._run_edge_step(
                    (parent_parent_node, parent_node), g
                )

        logger.debug(f"Running step of full dag edge: {edge}")
        g.get_edge_data(*edge)["step"].run()


class RunnableJobSequenceWithReport(RunnableJobWithReport):
    def __init__(self, job_dir_path, criterion, top_n_jobs_to_keep, steps):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
        )

        # validate
        for step in steps:
            validate_variable_type(step, RunnableJob)

        #
        self.steps = steps

    def run(self):
        # run steps in sequence
        for step in self.steps:
            step.run()

        # write results CSV
        self.write_results_csv()

        # write report PDF
        self.write_report_pdf()

    def write_results_csv(self):
        # TODO: write / overwrite special column identifying ordinal position in sequence
        df = pd.DataFrame()
        for step in self.steps:
            df = pd.concat([df, step.load_results_csv()], ignore_index=True)
        df.to_csv()


class Pipeline(RunnableJobSequenceWithReport):
    def __init__(self, job_dir_path, criterion, top_n_jobs_to_keep, steps):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
            steps=steps,
        )


class RepeatableRunnableSequence(RunnableJobSequenceWithReport):
    def __init__(
        self,
        job_dir_path,
        criterion,
        top_n_jobs_to_keep,
        steps,
        n_repetitions=0,
        max_iterations_with_no_improvement=sys.maxsize,
    ):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
            steps=steps,
        )

        #
        self.n_repetitions = n_repetitions
        self.max_iterations_with_no_improvement = max_iterations_with_no_improvement

    def run(self):
        # TODO
        pass


class RunnableSequenceIteration(RunnableJobSequenceWithReport):
    def __init__(self, job_dir_path, criterion, top_n_jobs_to_keep, steps):
        super().__init__(
            job_dir_path=job_dir_path,
            criterion=criterion,
            top_n_jobs_to_keep=top_n_jobs_to_keep,
            steps=steps,
        )
