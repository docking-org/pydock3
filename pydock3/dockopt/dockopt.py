from typing import Union, Iterable, List, Tuple
import itertools
import os
import shutil
import sys
from functools import wraps
from dataclasses import dataclass, fields, asdict
from copy import copy, deepcopy
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

import seaborn as sns
from joypy import joyplot


from pydock3.util import (
    Script,
    CleanExit,
    get_dataclass_as_dict,
    validate_variable_type,
    get_hexdigest_of_persistent_md5_hash_of_tuple,
)
from pydock3.config import (
    Parameter,
    flatten_and_parameter_cast_param_dict,
    get_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict,
)
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
    BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT,
    DOCK_FILE_IDENTIFIERS,
    WorkingDir,
    BlasterFile,
    DockFiles,
    BlasterStep,
)
from pydock3.dockopt.roc import ROC
from pydock3.jobs import ArrayDockingJob, DOCK3_EXECUTABLE_PATH
from pydock3.job_schedulers import SlurmJobScheduler, SGEJobScheduler
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH
from pydock3.retrodock.retrodock import log_job_submission_result, get_results_dataframe_from_actives_job_and_decoys_job_outdock_files, str_to_float, ROC_PLOT_FILE_NAME
from pydock3.blastermaster.util import DEFAULT_FILES_DIR_PATH
from pydock3.dockopt.results import BEST_RETRODOCK_JOBS_DIR_NAME, RESULTS_CSV_FILE_NAME, ResultsManager, DockoptStepResultsManager, DockoptStepSequenceIterationResultsManager, DockoptStepSequenceResultsManager
from pydock3.dockopt.reporter import Reporter, RETRODOCK_JOB_ID_COLUMN_NAME
from pydock3.dockopt.criterion import EnrichmentScore, Criterion
from pydock3.dockopt.pipeline import PipelineComponent, PipelineComponentSequence, PipelineComponentSequenceIteration
from pydock3.dockopt.parameters import DockoptComponentParametersManager
from pydock3.dockopt.docking_configuration import DockFileNodesTuple, DockingConfiguration
from pydock3.dockopt.dock_files_modification.matching_spheres_perturbation import MatchingSpheresPerturbationStep

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
SCHEDULER_NAME_TO_CLASS_DICT = {
    "sge": SGEJobScheduler,
    "slurm": SlurmJobScheduler,
}

#
CRITERION_CLASS_DICT = {"enrichment_score": EnrichmentScore}


@dataclass
class DockoptPipelineComponentRunFuncArgSet:
    scheduler: str
    actives_tgz_file_path: Union[None, str] = None
    decoys_tgz_file_path: Union[None, str] = None
    temp_storage_path: Union[None, str] = None
    retrodock_job_max_reattempts: int = 0
    retrodock_job_timeout_minutes: Union[None, int] = None
    max_scheduler_jobs_running_at_a_time: Union[None, int] = None
    export_decoy_poses: bool = False


class Dockopt(Script):
    JOB_DIR_NAME = "dockopt_job"
    CONFIG_FILE_NAME = "dockopt_config.yaml"
    ACTIVES_TGZ_FILE_NAME = "actives.tgz"
    DECOYS_TGZ_FILE_NAME = "decoys.tgz"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(
        os.path.dirname(DOCKOPT_INIT_FILE_PATH), "default_dockopt_config.yaml"
    )

    def __init__(self):
        super().__init__()

        #
        self.job_dir = None  # assigned in .init()

    @staticmethod
    def handle_run_func(run_func):
        @wraps(run_func)
        def wrapper(self, *args, **kwargs):
            with CleanExit():
                logger.info(f"Running {self.__class__.__name__}")
                run_func(self, *args, **kwargs)

        return wrapper

    def init(
        self,
        job_dir_path: str = JOB_DIR_NAME,
        overwrite: bool = False
    ) -> None:
        # create job dir
        self.job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create working dir & copy in blaster files
        blaster_file_names = list(BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.values())
        user_provided_blaster_file_paths = [
            os.path.abspath(f) for f in blaster_file_names if os.path.isfile(f)
        ]
        files_to_copy_str = "\n\t".join(user_provided_blaster_file_paths)
        if user_provided_blaster_file_paths:
            logger.info(
                f"Copying the following files from current directory into job working directory:\n\t{files_to_copy_str}"
            )
            for blaster_file_path in user_provided_blaster_file_paths:
                self.job_dir.copy_in_file(blaster_file_path)
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
        scheduler: str,
        job_dir_path: str = ".",
        config_file_path: str = None,
        actives_tgz_file_path: str = None,
        decoys_tgz_file_path: str = None,
        retrodock_job_max_reattempts: int = 0,
        retrodock_job_timeout_minutes: str = None,
        max_scheduler_jobs_running_at_a_time: str = None,  # TODO
        export_decoy_poses: bool = False,  # TODO
    ) -> None:
        # validate args
        if config_file_path is None:
            config_file_path = os.path.join(job_dir_path, self.CONFIG_FILE_NAME)
        if actives_tgz_file_path is None:
            actives_tgz_file_path = os.path.join(job_dir_path, self.ACTIVES_TGZ_FILE_NAME)
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
            temp_storage_path = os.environ["TMPDIR"]
        except KeyError:
            logger.error(
                "The following environmental variables are required to submit retrodock jobs: TMPDIR"
            )
            return

        #
        component_run_func_arg_set = DockoptPipelineComponentRunFuncArgSet(
            scheduler=scheduler,
            actives_tgz_file_path=actives_tgz_file_path,
            decoys_tgz_file_path=decoys_tgz_file_path,
            temp_storage_path=temp_storage_path,
            retrodock_job_max_reattempts=retrodock_job_max_reattempts,
            retrodock_job_timeout_minutes=retrodock_job_timeout_minutes,
            max_scheduler_jobs_running_at_a_time=max_scheduler_jobs_running_at_a_time,
            export_decoy_poses=export_decoy_poses,
        )

        #
        logger.info("Loading config file...")
        config = DockoptParametersConfiguration(config_file_path)
        logger.info("done.")

        #
        config_params_str = "\n".join(
            [
                f"{param_name}: {param.value}"
                for param_name, param in flatten_and_parameter_cast_param_dict(
                    config.param_dict
                ).items()
            ]
        )
        logger.debug(f"Parameters:\n{config_params_str}")

        #
        proper_blaster_file_names = list(BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.values())
        blaster_files_to_copy_in = [
            os.path.abspath(f) for f in proper_blaster_file_names if os.path.isfile(f)
        ]

        #
        pipeline = DockoptStepSequenceIteration(
            **config.param_dict["pipeline"],
            component_id=File.get_file_name_of_file(job_dir_path),
            dir_path=job_dir_path,
            blaster_files_to_copy_in=blaster_files_to_copy_in,
        )
        pipeline.run(component_run_func_arg_set=component_run_func_arg_set)


class DockoptStep(PipelineComponent):
    WORKING_DIR_NAME = "working"
    RETRODOCK_JOBS_DIR_NAME = "retrodock_jobs"

    def __init__(
            self,
            component_id: str,
            dir_path: str,
            criterion: str,
            top_n: int,
            parameters: Iterable[dict],
            dock_files_to_use_from_previous_component: dict,
            blaster_files_to_copy_in: Iterable[BlasterFile],
            last_component_completed: Union[PipelineComponent, None] = None,
    ):
        super().__init__(
            component_id=component_id,
            dir_path=dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=DockoptStepResultsManager(RESULTS_CSV_FILE_NAME),
        )

        #
        blaster_file_names = list(BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.values())
        backup_blaster_file_paths = [
            os.path.join(DEFAULT_FILES_DIR_PATH, blaster_file_name)
            for blaster_file_name in blaster_file_names
        ]
        new_file_names = [
            f"{File.get_file_name_of_file(file_path)}_1"  # TODO
            for file_path in blaster_files_to_copy_in
        ]  # all nodes in graph will be numerically indexed, including input files
        new_backup_file_names = [
            f"{File.get_file_name_of_file(file_path)}_1"  # TODO
            for file_path in backup_blaster_file_paths
        ]  # ^
        # TODO: remove need to copy these blaster files in for every step
        self.working_dir = WorkingDir(
            path=os.path.join(self.dir.path, self.WORKING_DIR_NAME),
            create=True,
            reset=False,
            files_to_copy_in=blaster_files_to_copy_in,
            new_file_names=new_file_names,
            backup_files_to_copy_in=backup_blaster_file_paths,
            new_backup_file_names=new_backup_file_names,
        )
        self.retrodock_jobs_dir = Dir(
            path=os.path.join(self.dir.path, self.RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=False,
        )

        #
        self.best_retrodock_jobs_dir = Dir(
            path=os.path.join(self.dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.actives_tgz_file = None  # set at beginning of .run()  # TODO: make non-class var
        self.decoys_tgz_file = None  # set at beginning of .run()  # TODO: make non-class var
        
        #
        if isinstance(parameters["dock_executable_path"], list):
            dock_executable_paths = []
            for dock_executable_path in parameters["dock_executable_path"].value:
                if dock_executable_path is None:
                    dock_executable_paths.append(DOCK3_EXECUTABLE_PATH)
                else:
                    dock_executable_paths.append(dock_executable_path)
        else:
            if parameters["dock_executable_path"] is None:
                dock_executable_paths = [DOCK3_EXECUTABLE_PATH]
            else:
                dock_executable_paths = [parameters["dock_executable_path"]]

        #
        dock_files_generation_flat_param_dicts = get_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["dock_files_generation"])
        dock_files_modification_flat_param_dicts = get_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["dock_files_modification"])
        indock_file_generation_flat_param_dicts = get_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["indock_file_generation"])

        # TODO: same thing happens in DockingConfiguration(). Maybe abstract this into another object?
        param_dict_hashes = []
        for p_dict in dock_files_generation_flat_param_dicts:
            p_dict_items_interleaved_sorted_by_key_tuple = tuple(
                itertools.chain.from_iterable(
                    sorted(list(zip(*list(zip(*p_dict.items())))), key=lambda x: x[0])
                )
            )
            param_dict_hashes.append(
                get_hexdigest_of_persistent_md5_hash_of_tuple(
                    p_dict_items_interleaved_sorted_by_key_tuple
                )
            )
        dock_files_generation_flat_param_dicts = [
            x
            for x, y in sorted(
                zip(dock_files_generation_flat_param_dicts, param_dict_hashes),
                key=lambda pair: pair[1],
            )
        ]

        #
        graph = nx.DiGraph()

        # for recording combinations of dock files
        partial_dock_file_identifier_to_node_id_dicts = []
        if last_component_completed is not None:
            for row_index, row in last_component_completed.load_results_dataframe().head(last_component_completed.top_n).iterrows():
                dock_file_identifier_to_node_id_dict = {}
                for dock_file_identifier, should_be_used in dock_files_to_use_from_previous_component.items():
                    if should_be_used:
                        #
                        dock_file_node_id = row.to_dict()[f"dockfiles.{dock_file_identifier}.node_id"]
                        dock_file_identifier_to_node_id_dict[dock_file_identifier] = dock_file_node_id
                        dock_file_lineage_subgraph = self._get_dock_file_lineage_subgraph(
                            graph=last_component_completed.graph,
                            dock_file_node_id=dock_file_node_id,
                        )
                        graph = nx.compose(graph, dock_file_lineage_subgraph)
                partial_dock_file_identifier_to_node_id_dicts.append(dock_file_identifier_to_node_id_dict)

        #
        dock_files_generation_dict_and_nodes_tuples = []
        dock_file_identifier_counter_dict = collections.defaultdict(int)
        dock_file_node_id_to_numerical_suffix_dict = {}
        blaster_files = BlasterFiles(working_dir=self.working_dir)
        for dock_files_generation_flat_param_dict in dock_files_generation_flat_param_dicts:
            # get config for get_blaster_steps
            # each value in dict must be an instance of Parameter
            steps = get_blaster_steps(
                blaster_files=blaster_files,
                flat_param_dict=dock_files_generation_flat_param_dict,
                working_dir=self.working_dir,
            )

            # form subgraph for this dock_files_generation_param_dict from the blaster steps it defines
            subgraph = self._get_graph_from_all_steps_in_order(self.component_id, steps)

            #
            dock_file_identifier_to_node_id_dict = {}
            for dock_file_identifier, should_be_used in dock_files_to_use_from_previous_component.items():
                step_hash_to_edges_dict = collections.defaultdict(list)
                step_hash_to_step_class_instance_dict = {}
                if not should_be_used:  # need to create during this dockopt step, so add to graph
                    #
                    dock_file_node_id = self._get_blaster_file_node_with_blaster_file_identifier(dock_file_identifier, subgraph)
                    dock_file_identifier_to_node_id_dict[dock_file_identifier] = dock_file_node_id
                    dock_file_lineage_subgraph = self._get_dock_file_lineage_subgraph(
                        graph=subgraph,
                        dock_file_node_id=dock_file_node_id,
                    ).copy()

                    #
                    new_dock_file_lineage_subgraph = deepcopy(dock_file_lineage_subgraph)
                    for node_id in self._get_blaster_file_nodes(dock_file_lineage_subgraph):
                        if node_id not in dock_file_node_id_to_numerical_suffix_dict:
                            blaster_file_identifier = dock_file_lineage_subgraph.nodes[node_id]['blaster_file'].identifier
                            dock_file_node_id_to_numerical_suffix_dict[node_id] = dock_file_identifier_counter_dict[blaster_file_identifier] + 1
                            dock_file_identifier_counter_dict[blaster_file_identifier] += 1
                        new_blaster_file = deepcopy(dock_file_lineage_subgraph.nodes[node_id]['blaster_file'])
                        new_blaster_file.path = f"{new_blaster_file.path}_{dock_file_node_id_to_numerical_suffix_dict[node_id]}"
                        new_dock_file_lineage_subgraph.nodes[node_id]['blaster_file'] = new_blaster_file
                    dock_file_lineage_subgraph = new_dock_file_lineage_subgraph

                    #
                    for u, v, data in dock_file_lineage_subgraph.edges(data=True):
                        step_hash_to_edges_dict[data["step_hash"]].append((u, v))
                        step_hash_to_step_class_instance_dict[data["step_hash"]] = data["step_instance"]

                    #
                    for step_hash, edges in step_hash_to_edges_dict.items():
                        #
                        kwargs = {"working_dir": self.working_dir}
                        for (parent_node, child_node) in edges:
                            edge_data_dict = dock_file_lineage_subgraph.get_edge_data(parent_node, child_node)
                            parent_node_data_dict = dock_file_lineage_subgraph.nodes[parent_node]
                            child_node_data_dict = dock_file_lineage_subgraph.nodes[child_node]
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
                        step_class = dock_file_lineage_subgraph.get_edge_data(*edges[0])["step_class"]  # first edge is fine since all edges have same step class
                        step_hash_to_step_class_instance_dict[step_hash] = step_class(**kwargs)

                        #
                        for parent_node, child_node in edges:
                            dock_file_lineage_subgraph.get_edge_data(parent_node, child_node)[
                                "step_instance"
                            ] = step_hash_to_step_class_instance_dict[step_hash]

                    #
                    graph = nx.compose(graph, new_dock_file_lineage_subgraph)

            #
            if partial_dock_file_identifier_to_node_id_dicts:
                for partial_dock_file_identifier_to_node_id_dict in partial_dock_file_identifier_to_node_id_dicts:
                    dock_file_nodes_tuple = DockFileNodesTuple(**{**partial_dock_file_identifier_to_node_id_dict, **dock_file_identifier_to_node_id_dict})  # complement + complement = complete
                    dock_files_generation_dict_and_nodes_tuples.append((dock_files_generation_flat_param_dict, dock_file_nodes_tuple))
            else:
                dock_file_nodes_tuple = DockFileNodesTuple(**dock_file_identifier_to_node_id_dict)  # use dock file nodes from this step only
                dock_files_generation_dict_and_nodes_tuples.append((dock_files_generation_flat_param_dict, dock_file_nodes_tuple))

        # matching spheres perturbation
        if dock_files_modification_flat_param_dicts[0][  # this should always hold
            "matching_spheres_perturbation.use"
        ].value:
            #
            unique_matching_spheres_file_nodes = list(set([x[1].matching_spheres_file for x in dock_files_generation_dict_and_nodes_tuples]))

            #
            perturbed_node_to_modification_dicts_index_dict = {}
            matching_spheres_node_to_perturbed_nodes_dict = collections.defaultdict(list)
            for dock_files_generation_dict, dock_file_nodes_tuple in dock_files_generation_dict_and_nodes_tuples:
                #
                for modification_dicts_index, dock_files_modification_flat_param_dict in enumerate(dock_files_modification_flat_param_dicts):
                    matching_spheres_blaster_file = graph.nodes[dock_file_nodes_tuple.matching_spheres_file]['blaster_file']
                    for i in range(
                        int(
                            dock_files_modification_flat_param_dict[
                                "matching_spheres_perturbation.num_samples_per_matching_spheres_file"
                            ].value
                        )
                    ):
                        max_deviation_angstroms = int(
                            dock_files_modification_flat_param_dict[
                                "matching_spheres_perturbation.max_deviation_angstroms"
                            ].value
                        )
                        max_deviation_angstroms_parameter = Parameter(
                            "matching_spheres_perturbation.max_deviation_angstroms",
                            max_deviation_angstroms,
                        )
                        perturbed_matching_spheres_file_path = os.path.join(
                            self.working_dir.path,
                            f"{matching_spheres_blaster_file.name}_{i+1}"  # file name will grow with number of perturbations
                        )
                        perturbed_matching_spheres_file = BlasterFile(perturbed_matching_spheres_file_path, identifier="matching_spheres_file")
                        step = MatchingSpheresPerturbationStep(
                            self.working_dir,
                            matching_spheres_infile=matching_spheres_blaster_file,
                            perturbed_matching_spheres_outfile=perturbed_matching_spheres_file,
                            max_deviation_angstroms_parameter=max_deviation_angstroms_parameter,
                        )

                        # get step hash from infile hashes, step dir, parameters, and outfiles
                        step_hash = DockoptStep._get_step_hash(component_id, step)

                        # add outfile node
                        outfile, = list(step.outfiles._asdict().values())
                        outfile_hash = DockoptStep._get_outfile_hash(component_id, outfile, step_hash)
                        graph.add_node(
                            outfile_hash,
                            blaster_file=deepcopy(outfile.original_file_in_working_dir),
                        )

                        # add parameter node
                        parameter, = list(step.parameters._asdict().values())
                        graph.add_node(parameter.hexdigest_of_persistent_md5_hash, parameter=deepcopy(parameter))

                        # connect each infile node to outfile node
                        infile_step_var_name, = list(step.infiles._asdict().keys())
                        outfile_step_var_name, = list(step.outfiles._asdict().keys())
                        parameter_step_var_name, = list(step.parameters._asdict().keys())
                        graph.add_edge(
                            dock_file_nodes_tuple.matching_spheres_file,
                            outfile_hash,
                            step_class=step.__class__,
                            original_step_dir_name=step.step_dir.name,
                            step_instance=deepcopy(step),
                            step_hash=step_hash,
                            parent_node_step_var_name=infile_step_var_name,
                            child_node_step_var_name=outfile_step_var_name,
                        )

                        # connect each parameter node to outfile node
                        graph.add_edge(
                            parameter.hexdigest_of_persistent_md5_hash,
                            outfile_hash,
                            step_class=step.__class__,
                            original_step_dir_name=step.step_dir.name,
                            step_instance=deepcopy(step),
                            step_hash=step_hash,
                            parent_node_step_var_name=parameter_step_var_name,
                            child_node_step_var_name=outfile_step_var_name,
                        )

                        #
                        matching_spheres_node_to_perturbed_nodes_dict[dock_file_nodes_tuple.matching_spheres_file].append(outfile_hash)
                        perturbed_node_to_modification_dicts_index_dict[outfile_hash] = modification_dicts_index

            #
            partial_docking_configuration_tuples = []
            for dock_files_generation_dict_and_nodes_tuple in dock_files_generation_dict_and_nodes_tuples:
                matching_spheres_file_node_id = dock_files_generation_dict_and_nodes_tuple[1].matching_spheres_file
                for perturbed_file_node_id in matching_spheres_node_to_perturbed_nodes_dict[matching_spheres_file_node_id]:
                    #
                    dock_files_generation_dict, dock_file_nodes_tuple = deepcopy(dock_files_generation_dict_and_nodes_tuple)
                    kwargs = {**{dock_file_identifier: getattr(dock_file_nodes_tuple, dock_file_identifier) for dock_file_identifier in DOCK_FILE_IDENTIFIERS}, **{'matching_spheres_file': perturbed_file_node_id}}
                    new_dock_file_nodes_tuple = DockFileNodesTuple(**kwargs)
                    dock_files_modification_dict = dock_files_modification_flat_param_dicts[perturbed_node_to_modification_dicts_index_dict[perturbed_file_node_id]]
                    partial_docking_configuration_tuple = (
                        dock_files_generation_dict,
                        dock_files_modification_dict,
                        new_dock_file_nodes_tuple,
                    )
                    partial_docking_configuration_tuples.append(partial_docking_configuration_tuple)
        else:
            partial_docking_configuration_tuples = []
            for dock_files_generation_dict_and_nodes_tuple in dock_files_generation_dict_and_nodes_tuples:
                dock_files_generation_dict, dock_file_nodes_tuple = deepcopy(dock_files_generation_dict_and_nodes_tuple)
                dock_files_modification_dict = dock_files_modification_flat_param_dicts[0]  # this should always hold (see above, same comment)
                partial_docking_configuration_tuple = (
                    dock_files_generation_dict,
                    dock_files_modification_dict,
                    dock_file_nodes_tuple,
                )
                partial_docking_configuration_tuples.append(partial_docking_configuration_tuple)
        #
        docking_configurations = []
        for dock_files_generation_flat_param_dict, dock_files_modification_flat_param_dict, dock_file_nodes_tuple in partial_docking_configuration_tuples:
            for dock_executable_path in dock_executable_paths:
                for indock_file_generation_flat_param_dict in indock_file_generation_flat_param_dicts:
                    configuration_num = len(docking_configurations) + 1
                    docking_configuration = DockingConfiguration(
                        component_id=self.component_id,
                        configuration_num=configuration_num,
                        dock_executable_path=dock_executable_path,
                        dock_files_generation_flat_param_dict=dock_files_generation_flat_param_dict,
                        dock_files_modification_flat_param_dict=dock_files_modification_flat_param_dict,
                        indock_file_generation_flat_param_dict=indock_file_generation_flat_param_dict,
                        dock_file_nodes_tuple=dock_file_nodes_tuple,
                        dock_files=DockFiles(**{dock_file_identifier: graph.nodes[node_id]['blaster_file'] for dock_file_identifier, node_id in dock_file_nodes_tuple._asdict().items()}),
                        indock_file=IndockFile(path=os.path.join(self.working_dir.path, f"{INDOCK_FILE_NAME}_{configuration_num}")),
                    )
                    docking_configurations.append(docking_configuration)

        # validate that there are no cycles (i.e. that it is a directed acyclic graph)
        if not nx.is_directed_acyclic_graph(graph):
            raise Exception("Cycle found in graph!")

        #
        self.graph = graph
        logger.debug(
            f"Graph initialized with:\n\tNodes: {self.graph.nodes}\n\tEdges: {self.graph.edges}"
        )

        #
        self.docking_configurations = docking_configurations

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.core.frame.DataFrame:
        if component_run_func_arg_set.actives_tgz_file_path is not None:
            self.actives_tgz_file = File(path=component_run_func_arg_set.actives_tgz_file_path)
        else:
            self.actives_tgz_file = None
        if component_run_func_arg_set.decoys_tgz_file_path is not None:
            self.decoys_tgz_file = File(path=component_run_func_arg_set.decoys_tgz_file_path)
        else:
            self.decoys_tgz_file = None

        # run necessary steps to get all dock files
        logger.info("Generating dock files & INDOCK for all docking configurations...")
        for docking_configuration in self.docking_configurations:
            # make dock files
            for dock_file_identifier in DOCK_FILE_IDENTIFIERS:
                self._run_unrun_steps_needed_to_create_this_blaster_file_node(
                    getattr(docking_configuration.dock_file_nodes_tuple, dock_file_identifier), self.graph
                )

            # make indock file now that dock files exist
            docking_configuration.indock_file.write(docking_configuration.dock_files, docking_configuration.full_flat_parameters_dict)
        logger.info("done.")

        #
        job_hash = get_hexdigest_of_persistent_md5_hash_of_tuple(tuple(sorted([docking_configuration.hexdigest_of_persistent_md5_hash for docking_configuration in self.docking_configurations])))
        with open(os.path.join(self.dir.path, "job_hash_hexdigest.md5"), "w") as f:
            f.write(f"{job_hash}\n")

        #
        array_job_docking_configurations_file_path = os.path.join(self.dir.path, "array_job_docking_configurations.txt")
        with open(array_job_docking_configurations_file_path, 'w') as f:
            for docking_configuration in self.docking_configurations:
                dockfile_paths_str = " ".join([self.graph.nodes[getattr(docking_configuration.dock_file_nodes_tuple, dock_file_identifier)]['blaster_file'].path for dock_file_identifier in DOCK_FILE_IDENTIFIERS])
                indock_file_path_str = docking_configuration.indock_file.path
                f.write(f"{docking_configuration.configuration_num} {indock_file_path_str} {dockfile_paths_str} {docking_configuration.dock_executable_path}\n")

        #
        def get_actives_outdock_file_path_for_configuration_num(configuration_num):
            return os.path.join(self.retrodock_jobs_dir.path, 'actives', str(configuration_num), 'OUTDOCK.0')

        def get_decoys_outdock_file_path_for_configuration_num(configuration_num):
            return os.path.join(self.retrodock_jobs_dir.path, 'decoys', str(configuration_num), 'OUTDOCK.0')

        # submit retrodock jobs (one for actives, one for decoys)
        array_jobs = []
        for sub_dir_name, should_export_mol2, input_molecules_tgz_file_path in [
            ('actives', True, component_run_func_arg_set.actives_tgz_file_path),
            ('decoys', False, component_run_func_arg_set.decoys_tgz_file_path),
        ]:
            array_job = ArrayDockingJob(
                name=f"dockopt_job_{job_hash}_{sub_dir_name}",
                job_dir=Dir(os.path.join(self.retrodock_jobs_dir.path, sub_dir_name)),
                input_molecules_tgz_file_path=input_molecules_tgz_file_path,
                job_scheduler=component_run_func_arg_set.scheduler,
                temp_storage_path=component_run_func_arg_set.temp_storage_path,
                array_job_docking_configurations_file_path=array_job_docking_configurations_file_path,
                max_reattempts=component_run_func_arg_set.retrodock_job_max_reattempts,
                export_mol2=should_export_mol2,
            )
            sub_result, procs = array_job.submit_all_tasks(
                job_timeout_minutes=component_run_func_arg_set.retrodock_job_timeout_minutes,
                skip_if_complete=True,
            )
            array_jobs.append(array_job)
            log_job_submission_result(array_job, sub_result, procs)

        # make a queue of tuples containing job-relevant data for processing
        docking_configurations_processing_queue = deepcopy(self.docking_configurations)

        # process results of docking jobs
        logger.info(
            f"Awaiting / processing retrodock job results ({len(docking_configurations_processing_queue)} tasks in total)"
        )
        data_dicts = []
        configuration_num_to_num_reattempts_dict = collections.defaultdict(int)
        while len(docking_configurations_processing_queue) > 0:
            #
            docking_configuration = docking_configurations_processing_queue.pop(0)

            if any([not array_job.task_is_complete(str(docking_configuration.configuration_num)) for array_job in array_jobs]):  # one or both OUTDOCK files do not exist yet
                time.sleep(
                    0.01
                )  # sleep for a bit
                if any([(not array_job.task_is_complete(str(docking_configuration.configuration_num))) and (not array_job.is_running) for array_job in array_jobs]):
                    # task must have timed out / failed for one or both jobs
                    logger.warning(
                        f"Failure / time out witnessed for task {docking_configuration.configuration_num}"
                    )
                    if configuration_num_to_num_reattempts_dict[docking_configuration.configuration_num] > component_run_func_arg_set.retrodock_job_max_reattempts:
                        logger.warning(
                            f"Max reattempts exhausted for task {docking_configuration.configuration_num}"
                        )
                        continue  # move on to next in queue without re-attempting failed task

                    for array_job in array_jobs:
                        if not array_job.task_is_complete(str(docking_configuration.configuration_num)):
                            array_job.submit_task(
                                str(docking_configuration.configuration_num),
                                job_timeout_minutes=component_run_func_arg_set.retrodock_job_timeout_minutes,
                                skip_if_complete=False,
                            )  # re-attempt relevant job(s) for incomplete task
                    configuration_num_to_num_reattempts_dict[docking_configuration.configuration_num] += 1

                docking_configurations_processing_queue.append(
                    docking_configuration
                )  # move to back of queue
                continue  # move on to next in queue

            # load outdock files and get dataframe
            try:
                # get dataframe of actives job results and decoys job results combined
                df = (
                    get_results_dataframe_from_actives_job_and_decoys_job_outdock_files(
                        get_actives_outdock_file_path_for_configuration_num(docking_configuration.configuration_num), get_decoys_outdock_file_path_for_configuration_num(docking_configuration.configuration_num)
                    )
                )
            except Exception as e:  # if outdock files failed to be parsed then re-attempt task
                logger.warning(f"Failed to parse outdock file(s) due to error: {e}")
                if configuration_num_to_num_reattempts_dict[docking_configuration.configuration_num] > component_run_func_arg_set.retrodock_job_max_reattempts:
                    logger.warning(
                        f"Max reattempts exhausted for task {docking_configuration.configuration_num}"
                    )
                    continue  # move on to next in queue without re-attempting failed task

                for array_job in array_jobs:
                    array_job.submit_task(
                        str(docking_configuration.configuration_num),
                        job_timeout_minutes=component_run_func_arg_set.retrodock_job_timeout_minutes,
                        skip_if_complete=False,
                    )  # re-attempt both jobs
                configuration_num_to_num_reattempts_dict[docking_configuration.configuration_num] += 1
                docking_configurations_processing_queue.append(
                    docking_configuration
                )  # move to back of queue
                continue  # move on to next in queue

            #
            logger.info(
                f"Task {docking_configuration.configuration_num} completed. Successfully loaded both OUTDOCK files."
            )

            # sort dataframe by total energy score
            df["total_energy"] = df["total_energy"].astype(float)
            df = df.sort_values(
                by=["total_energy", "is_active"], na_position="last", ignore_index=True
            )  # sorting secondarily by 'is_active' (0 or 1) ensures that decoys are ranked before actives in case they have the same exact score (pessimistic approach)
            df = df.drop_duplicates(
                subset=["db2_file_path"], keep="first", ignore_index=True
            )

            # make data dict for this job (will be used to make dataframe for results of all jobs)
            data_dict = {f"parameters.{key}": value for key, value in docking_configuration.full_flat_parameters_dict.items()}
            data_dict[RETRODOCK_JOB_ID_COLUMN_NAME] = str(docking_configuration.configuration_num)

            #
            for dock_file_identifier in DOCK_FILE_IDENTIFIERS:
                dock_file_node_id = getattr(docking_configuration.dock_file_nodes_tuple, dock_file_identifier)
                dock_file = self.graph.nodes[dock_file_node_id]['blaster_file']
                data_dict[f"dockfiles.{dock_file_identifier}.name"] = dock_file.name
                data_dict[f"dockfiles.{dock_file_identifier}.node_id"] = dock_file_node_id

            #
            data_dict["dockfiles.indock_file.name"] = docking_configuration.indock_file.name

            # get ROC and calculate enrichment score of this job's docking set-up
            if isinstance(self.criterion, EnrichmentScore):
                logger.debug("Calculating ROC and enrichment score...")
                booleans = df["is_active"]
                data_dict[self.criterion.name] = self.criterion.calculate(booleans)
                logger.debug("done.")

            #
            data_dict['pipeline_component_id'] = self.component_id

            # save data_dict for this job
            data_dicts.append(data_dict)

        # write jobs completion status
        num_tasks_successful = len(data_dicts)
        logger.info(
            f"Finished {num_tasks_successful} out of {len(self.docking_configurations)} tasks."
        )
        if num_tasks_successful == 0:
            raise Exception(
                "All tasks failed. Something is wrong. Please check logs."
            )

        # make dataframe of optimization job results
        df = pd.DataFrame(data=data_dicts)

        return df

    @staticmethod
    def _get_dock_file_lineage_subgraph(graph, dock_file_node_id):
        """Gets the subgraph representing the steps necessary to produce the desired dock file"""

        node_ids = [dock_file_node_id] + list(nx.ancestors(graph, dock_file_node_id))
        step_hashes = []
        for u, v, data in graph.edges(data=True):
            if u in node_ids and v in node_ids:
                step_hashes.append(data['step_hash'])
        for u, v, data in graph.edges(data=True):
            if data['step_hash'] in step_hashes:
                node_ids.append(u)
                node_ids.append(v)
        node_ids = list(set(node_ids))

        return graph.subgraph(node_ids)

    @staticmethod
    def _get_infile_hash(component_id, infile):
        return get_hexdigest_of_persistent_md5_hash_of_tuple((component_id, infile.original_file_in_working_dir.name))

    @staticmethod
    def _get_outfile_hash(component_id, outfile, step_hash):
        return get_hexdigest_of_persistent_md5_hash_of_tuple((component_id, outfile.original_file_in_working_dir.name, step_hash))

    @staticmethod
    def _get_step_hash(component_id, step):
        #
        infiles_dict_items_list = sorted(step.infiles._asdict().items())
        outfiles_dict_items_list = sorted(step.outfiles._asdict().items())
        parameters_dict_items_list = sorted(step.parameters._asdict().items())

        #
        infile_hash_tuples = []
        for infile_step_var_name, infile in infiles_dict_items_list:
            infile_hash = DockoptStep._get_infile_hash(component_id, infile)
            infile_hash_tuples.append(
                (
                    infile_step_var_name,
                    infile_hash,
                )
            )
        return get_hexdigest_of_persistent_md5_hash_of_tuple(
            tuple(
                infile_hash_tuples
                + [step.__class__.__name__, step.step_dir.name]
                + [(parameter_step_var_name, parameter.hexdigest_of_persistent_md5_hash) for parameter_step_var_name, parameter in parameters_dict_items_list]
                + [(outfile_step_var_name, outfile.original_file_in_working_dir.name) for outfile_step_var_name, outfile in outfiles_dict_items_list]
            )
        )

    @staticmethod
    def _get_graph_from_all_steps_in_order(component_id: str, steps: List[BlasterStep]):
        #
        graph = nx.DiGraph()
        blaster_file_hash_dict = {}
        for step in steps:
            #
            infiles_dict_items_list = sorted(step.infiles._asdict().items())
            outfiles_dict_items_list = sorted(step.outfiles._asdict().items())
            parameters_dict_items_list = sorted(step.parameters._asdict().items())
    
            # get step hash from infile hashes, step dir, parameters, and outfiles
            step_hash = DockoptStep._get_step_hash(component_id, step)
    
            # add infile nodes
            for infile in step.infiles:
                if DockoptStep._get_blaster_file_node_with_same_file_name(infile.original_file_in_working_dir.name, graph) is not None:
                    continue
                if (component_id, infile.original_file_in_working_dir.name) not in blaster_file_hash_dict:
                    blaster_file_hash_dict[(component_id, infile.original_file_in_working_dir.name)] = DockoptStep._get_infile_hash(component_id, infile)
                graph.add_node(
                    blaster_file_hash_dict[(component_id, infile.original_file_in_working_dir.name)],
                    blaster_file=deepcopy(infile.original_file_in_working_dir),
                )
    
            # add outfile nodes
            for outfile_step_var_name, outfile in outfiles_dict_items_list:
                if DockoptStep._get_blaster_file_node_with_same_file_name(
                        outfile.original_file_in_working_dir.name, graph
                ):
                    raise Exception(
                        f"Attempting to add outfile to graph that already has said outfile as node: {outfile.original_file_in_working_dir.name}"
                    )
                if (component_id, outfile.original_file_in_working_dir.name) not in blaster_file_hash_dict:
                    blaster_file_hash_dict[(component_id, outfile.original_file_in_working_dir.name)] = DockoptStep._get_outfile_hash(component_id, outfile, step_hash)
                graph.add_node(
                    blaster_file_hash_dict[(component_id, outfile.original_file_in_working_dir.name)],
                    blaster_file=deepcopy(outfile.original_file_in_working_dir),
                )
    
            # add parameter nodes
            for parameter in step.parameters:
                graph.add_node(parameter.hexdigest_of_persistent_md5_hash, parameter=deepcopy(parameter))
    
            # connect each infile node to every outfile node
            for (infile_step_var_name, infile), (outfile_step_var_name, outfile) in itertools.product(
                infiles_dict_items_list, outfiles_dict_items_list
            ):
                graph.add_edge(
                    blaster_file_hash_dict[(component_id, infile.original_file_in_working_dir.name)],
                    blaster_file_hash_dict[(component_id, outfile.original_file_in_working_dir.name)],
                    step_class=step.__class__,
                    original_step_dir_name=step.step_dir.name,
                    step_instance=deepcopy(step),
                    step_hash=step_hash,
                    parent_node_step_var_name=infile_step_var_name,
                    child_node_step_var_name=outfile_step_var_name,
                )
    
            # connect each parameter node to every outfile nodes
            for (parameter_step_var_name, parameter), (outfile_step_var_name, outfile) in itertools.product(
                parameters_dict_items_list, outfiles_dict_items_list
            ):
                graph.add_edge(
                    parameter.hexdigest_of_persistent_md5_hash,
                    blaster_file_hash_dict[(component_id, outfile.original_file_in_working_dir.name)],
                    step_class=step.__class__,
                    original_step_dir_name=step.step_dir.name,
                    step_instance=deepcopy(step),  # this will be replaced with step instance with unique dir path
                    step_hash=step_hash,
                    parent_node_step_var_name=parameter_step_var_name,
                    child_node_step_var_name=outfile_step_var_name,
                )
    
        return graph

    @staticmethod
    def _get_blaster_file_nodes(g: nx.classes.digraph.DiGraph) -> str:
        return [node_id for node_id, node_data in g.nodes.items() if g.nodes[node_id].get("blaster_file")]
    
    @staticmethod
    def _get_blaster_file_node_with_blaster_file_identifier(
        blaster_file_identifier: str,
        g: nx.classes.digraph.DiGraph,
    ) -> str:
        blaster_file_node_ids = DockoptStep._get_blaster_file_nodes(g)
        if len(blaster_file_node_ids) == 0:
            return None
        matching_blaster_file_nodes = [node_id for node_id in blaster_file_node_ids if g.nodes[node_id]["blaster_file"].identifier == blaster_file_identifier]
        if len(matching_blaster_file_nodes) == 0:
            return None
        matching_blaster_file_node, = matching_blaster_file_nodes

        return matching_blaster_file_node

    @staticmethod
    def _get_blaster_file_node_with_same_file_name(
        file_name: str,
        g: nx.classes.digraph.DiGraph,
    ) -> BlasterFile:
        blaster_file_node_ids = DockoptStep._get_blaster_file_nodes(g)
        if len(blaster_file_node_ids) == 0:
            return None
        matching_blaster_file_nodes = [node_id for node_id in blaster_file_node_ids if file_name == g.nodes[node_id]["blaster_file"].name]
        if len(matching_blaster_file_nodes) == 0:
            return None
        matching_blaster_file_node, = matching_blaster_file_nodes

        return matching_blaster_file_node

    @staticmethod
    def _get_start_nodes(g: nx.classes.digraph.DiGraph) -> List[BlasterFile]:
        start_nodes = []
        for node_id, node_data in g.nodes.items():
            if g.in_degree(node_id) == 0:
                start_nodes.append(node_id)
        return start_nodes

    @staticmethod
    def _get_dock_file_nodes_tuple(g: nx.classes.digraph.DiGraph) -> List[str]:
        return DockFileNodesTuple(**{dock_file_identifier: DockoptStep._get_blaster_file_node_with_blaster_file_identifier(dock_file_identifier, g) for dock_file_identifier in DOCK_FILE_IDENTIFIERS})

    @staticmethod
    def _run_unrun_steps_needed_to_create_this_blaster_file_node(
        blaster_file_node: str,
        g: nx.classes.digraph.DiGraph,
    ):
        if g.nodes[blaster_file_node].get("blaster_file") is not None:
            blaster_file = g.nodes[blaster_file_node]['blaster_file']
            if not blaster_file.exists:
                for parent_node in g.predecessors(blaster_file_node):
                    DockoptStep._run_unrun_steps_needed_to_create_this_blaster_file_node(
                        parent_node, g
                    )
                a_parent_node = list(g.predecessors(blaster_file_node))[0]
                step_instance = g[a_parent_node][blaster_file_node]["step_instance"]
                if step_instance.is_done:
                    raise Exception(
                        f"blaster file {blaster_file.path} does not exist but step instance is_done=True"
                    )
                step_instance.run()


class DockoptStepSequenceIteration(PipelineComponentSequenceIteration):

    def __init__(
        self,
        component_id: str,
        dir_path: str,
        criterion: str,
        top_n: int,
        components: Iterable[dict],
        blaster_files_to_copy_in: Iterable[BlasterFile],
        last_component_completed: Union[PipelineComponent, None] = None,
        **kwargs  # TODO: make this unnecessary
    ):
        super().__init__(
            component_id=component_id,
            dir_path=dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=DockoptStepSequenceIterationResultsManager(RESULTS_CSV_FILE_NAME),
            components=components,
        )

        #
        self.blaster_files_to_copy_in = blaster_files_to_copy_in
        self.last_component_completed = last_component_completed

        #
        self.best_retrodock_jobs_dir = Dir(
            path=os.path.join(self.dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.core.frame.DataFrame:
        df = pd.DataFrame()
        best_criterion_value_witnessed = -float('inf')
        last_component_completed_in_sequence = self.last_component_completed
        for i, component_identifier_dict in enumerate(self.components):
            #
            component_num = i + 1

            #
            if "step" in component_identifier_dict:
                component_identifier = "step"
                component_class = DockoptStep
            elif "sequence" in component_identifier_dict:
                component_identifier = "sequence"
                component_class = DockoptStepSequence
            else:
                raise Exception(f"Dict must have one of 'step' or 'sequence' as keys. Witnessed: {component_identifier_dict}")

            #
            component_parameters_dict = deepcopy(component_identifier_dict[component_identifier])
            parameters_manager = DockoptComponentParametersManager(
                parameters_dict=component_parameters_dict,
                last_component_completed=last_component_completed_in_sequence,
            )

            #
            kwargs = {**parameters_manager.parameters_dict, **{
                'component_id': f"{self.component_id}.{component_num}",
                'dir_path': os.path.join(self.dir.path, str(component_num)),
                'blaster_files_to_copy_in': self.blaster_files_to_copy_in,  # TODO: is this necessary?
                'last_component_completed': last_component_completed_in_sequence,
            }}

            #
            component = component_class(**kwargs)

            #
            component.run(component_run_func_arg_set)

            #
            df_component = component.load_results_dataframe()
            df = pd.concat([df, df_component], ignore_index=True)

            #
            last_component_completed_in_sequence = component

        return df


class DockoptStepSequence(PipelineComponentSequence):

    def __init__(
        self,
        component_id: str,
        dir_path: str,
        criterion: str,
        top_n: int,
        components: Iterable[dict],
        num_repetitions: int,
        max_iterations_with_no_improvement: int,
        inter_iteration_criterion: str,
        inter_iteration_top_n: int,
        blaster_files_to_copy_in: Iterable[BlasterFile],
        last_component_completed: Union[PipelineComponent, None] = None,
        **kwargs  # TODO: make this unnecessary
    ):
        super().__init__(
            component_id=component_id,
            dir_path=dir_path,
            criterion=criterion,
            top_n=top_n,
            results_manager=DockoptStepSequenceResultsManager(RESULTS_CSV_FILE_NAME),
            components=components,
            num_repetitions=num_repetitions,
            max_iterations_with_no_improvement=max_iterations_with_no_improvement,
            inter_iteration_criterion=inter_iteration_criterion,
            inter_iteration_top_n=inter_iteration_top_n,
        )

        #
        self.blaster_files_to_copy_in = blaster_files_to_copy_in
        self.last_component_completed = last_component_completed

        #
        self.best_retrodock_jobs_dir = Dir(
            path=os.path.join(self.dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.core.frame.DataFrame:
        df = pd.DataFrame()
        best_criterion_value_witnessed = -float('inf')
        last_component_completed_in_sequence = self.last_component_completed
        num_iterations_left_with_no_improvement = self.max_iterations_with_no_improvement
        for i in range(self.num_repetitions+1):
            #
            iteration_num =  i + 1

            #
            if i == 0:
                component_criterion = self.last_component_completed.criterion.name
                component_top_n = self.last_component_completed.top_n
            else:
                component_criterion = self.inter_iteration_criterion
                component_top_n = self.inter_iteration_top_n

            #
            component = DockoptStepSequenceIteration(
                component_id=f"{self.component_id}.iter={iteration_num}",
                dir_path=os.path.join(self.dir.path, f"iter={iteration_num}"),
                criterion=component_criterion,
                top_n=component_top_n,
                components=self.components,
                blaster_files_to_copy_in=self.blaster_files_to_copy_in,
                last_component_completed=last_component_completed_in_sequence,
            )

            #
            component.run(component_run_func_arg_set)

            #
            df_component = component.load_results_dataframe()
            df = pd.concat([df, df_component], ignore_index=True)

            #
            last_component_completed_in_sequence = component

            #
            best_criterion_value_witnessed_this_iteration = df_component[component.criterion.name].max()
            if best_criterion_value_witnessed_this_iteration <= best_criterion_value_witnessed:
                if num_iterations_left_with_no_improvement == 0:
                    break
                else:
                    num_iterations_left_with_no_improvement -= 1
            else:
                best_criterion_value_witnessed = best_criterion_value_witnessed_this_iteration

        return df
