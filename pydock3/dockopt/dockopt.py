import uuid
from typing import Union, Iterable, List, Tuple, Callable, Any, TypeVar
from typing_extensions import ParamSpec
import itertools
import os
import shutil
import sys
from functools import wraps
from dataclasses import dataclass, fields, asdict, astuple
from copy import copy, deepcopy
import logging
import collections
import time
import random
from datetime import datetime, timedelta

import networkx as nx
import numpy as np
import pandas as pd
from dirhash import dirhash
import matplotlib.pyplot as plt

import seaborn as sns
from joypy import joyplot


from pydock3.util import (
    T,
    P,
    filter_kwargs_for_callable,
    Script,
    CleanExit,
    get_dataclass_as_dict,
    validate_variable_type,
    get_hexdigest_of_persistent_md5_hash_of_tuple,
)
from pydock3.dockopt.util import WORKING_DIR_NAME, RETRODOCK_JOBS_DIR_NAME, RESULTS_CSV_FILE_NAME, BEST_RETRODOCK_JOBS_DIR_NAME
from pydock3.config import (
    Parameter,
    flatten_and_parameter_cast_param_dict,
    get_sorted_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict,
)
from pydock3.blastermaster.blastermaster import BlasterFiles, get_blaster_steps
from pydock3.dockopt.config import DockoptParametersConfiguration
from pydock3.files import (
    INDOCK_FILE_NAME,
    Dir,
    File,
    IndockFile,
    OutdockFile,
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
from pydock3.jobs import ArrayDockingJob
from pydock3.job_schedulers import SlurmJobScheduler, SGEJobScheduler
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH
from pydock3.retrodock.retrodock import log_job_submission_result, get_results_dataframe_from_actives_job_and_decoys_job_outdock_files, str_to_float, ROC_PLOT_FILE_NAME
from pydock3.blastermaster.util import DEFAULT_FILES_DIR_PATH
from pydock3.dockopt.results import ResultsManager, DockoptStepResultsManager, DockoptStepSequenceIterationResultsManager, DockoptStepSequenceResultsManager
from pydock3.dockopt.reporter import Reporter
from pydock3.dockopt.criterion import EnrichmentScore, Criterion
from pydock3.dockopt.pipeline import PipelineComponent, PipelineComponentSequence, PipelineComponentSequenceIteration, Pipeline
from pydock3.dockopt.parameters import DockoptComponentParametersManager
from pydock3.dockopt.docking_configuration import DockingConfiguration, DockFileCoordinates, DockFileCoordinate, IndockFileCoordinate
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

#
MIN_SECONDS_BETWEEN_QUEUE_CHECKS = 2


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

    def __init__(self) -> None:
        super().__init__()

        #
        self.job_dir = None  # assigned in .init()

    @staticmethod
    def handle_run_func(run_func: Callable[P, T]) -> Callable[P, T]:
        @wraps(run_func)
        def wrapper(self, *args: P.args, **kwargs: P.kwargs):
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
        config_file_path: Union[str, None] = None,
        actives_tgz_file_path: Union[str, None] = None,
        decoys_tgz_file_path: Union[str, None] = None,
        retrodock_job_max_reattempts: int = 0,
        retrodock_job_timeout_minutes: Union[str, None] = None,
        max_scheduler_jobs_running_at_a_time: Union[str, None] = None,  # TODO
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
        pipeline = DockoptPipeline(
            **config.param_dict["pipeline"],
            pipeline_dir_path=job_dir_path,
            blaster_files_to_copy_in=blaster_files_to_copy_in,
        )
        pipeline.run(component_run_func_arg_set=component_run_func_arg_set)


class DockoptStep(PipelineComponent):
    def __init__(
            self,
            pipeline_dir_path: str,
            component_id: str,
            criterion: str,
            top_n: int,
            parameters: Iterable[dict],
            dock_files_to_use_from_previous_component: dict,
            blaster_files_to_copy_in: Iterable[BlasterFile],
            last_component_completed: Union[PipelineComponent, None] = None,
    ) -> None:
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
            component_id=component_id,
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
            path=os.path.join(self.component_dir.path, WORKING_DIR_NAME),
            create=True,
            reset=False,
            files_to_copy_in=blaster_files_to_copy_in,
            new_file_names=new_file_names,
            backup_files_to_copy_in=backup_blaster_file_paths,
            new_backup_file_names=new_backup_file_names,
        )
        self.retrodock_jobs_dir = Dir(
            path=os.path.join(self.component_dir.path, RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=False,
        )

        #
        self.best_retrodock_jobs_dir = Dir(
            path=os.path.join(self.component_dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.actives_tgz_file = None  # set at beginning of .run()  # TODO: make non-class var
        self.decoys_tgz_file = None  # set at beginning of .run()  # TODO: make non-class var

        #
        if isinstance(parameters["custom_dock_executable"], list):
            custom_dock_executables = [custom_dock_executable for custom_dock_executable in parameters["custom_dock_executable"]]
        else:
            custom_dock_executables = [parameters["custom_dock_executable"]]

        #
        sorted_dock_files_generation_flat_param_dicts = get_sorted_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["dock_files_generation"])
        sorted_dock_files_modification_flat_param_dicts = get_sorted_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["dock_files_modification"])
        sorted_indock_file_generation_flat_param_dicts = get_sorted_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(parameters["indock_file_generation"])

        #
        logger.debug(f"{len(sorted_dock_files_generation_flat_param_dicts)} dock file generation parametrizations:\n{sorted_dock_files_generation_flat_param_dicts}")
        logger.debug(f"{len(sorted_dock_files_modification_flat_param_dicts)} dock file modification parametrizations:\n{sorted_dock_files_modification_flat_param_dicts}")
        logger.debug(f"{len(sorted_indock_file_generation_flat_param_dicts)} indock file generation parametrizations:\n{sorted_indock_file_generation_flat_param_dicts}")

        #
        graph = nx.DiGraph()

        #
        last_component_docking_configurations = []
        if last_component_completed is not None and any(list(dock_files_to_use_from_previous_component.values())):
            logger.debug(f"Using the following dock files from previous component: {sorted([key for key, value in dock_files_to_use_from_previous_component.items() if value])}")
            for row_index, row in last_component_completed.load_results_dataframe().head(last_component_completed.top_n).iterrows():
                dc = DockingConfiguration.from_dict(row.to_dict())
                for dock_file_identifier, should_be_used in dock_files_to_use_from_previous_component.items():
                    if should_be_used:
                        #
                        dock_file_node_id = getattr(dc.dock_file_coordinates, dock_file_identifier).node_id
                        dock_file_lineage_subgraph = self._get_dock_file_lineage_subgraph(
                            graph=last_component_completed.graph,
                            dock_file_node_id=dock_file_node_id,
                        )
                        graph = nx.compose(graph, dock_file_lineage_subgraph)
                last_component_docking_configurations.append(dc)

        #
        dock_file_identifier_counter_dict = collections.defaultdict(int)
        blaster_file_node_id_to_numerical_suffix_dict = {}
        blaster_files = BlasterFiles(working_dir=self.working_dir)
        partial_dock_file_nodes_combination_dicts = []
        if any([not x for x in dock_files_to_use_from_previous_component.values()]):
            logger.debug(f"The following dock files will be generated during this step: {sorted([key for key, value in dock_files_to_use_from_previous_component.items() if not value])}")
            for dock_files_generation_flat_param_dict in sorted_dock_files_generation_flat_param_dicts:
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
                partial_dock_file_nodes_combination_dict = {}
                for dock_file_identifier, should_be_used in dock_files_to_use_from_previous_component.items():
                    step_hash_to_edges_dict = collections.defaultdict(list)
                    step_hash_to_step_class_instance_dict = {}
                    if not should_be_used:  # need to create during this dockopt step, so add to graph
                        #
                        dock_file_node_id = self._get_blaster_file_node_with_blaster_file_identifier(dock_file_identifier, subgraph)
                        partial_dock_file_nodes_combination_dict[dock_file_identifier] = dock_file_node_id
                        dock_file_lineage_subgraph = self._get_dock_file_lineage_subgraph(
                            graph=subgraph,
                            dock_file_node_id=dock_file_node_id,
                        )

                        #
                        new_dock_file_lineage_subgraph = deepcopy(dock_file_lineage_subgraph)
                        for node_id in self._get_blaster_file_nodes(dock_file_lineage_subgraph):
                            if node_id not in blaster_file_node_id_to_numerical_suffix_dict:
                                blaster_file_identifier = dock_file_lineage_subgraph.nodes[node_id]['blaster_file'].identifier
                                blaster_file_node_id_to_numerical_suffix_dict[node_id] = dock_file_identifier_counter_dict[blaster_file_identifier] + 1
                                dock_file_identifier_counter_dict[blaster_file_identifier] += 1
                            new_blaster_file = deepcopy(dock_file_lineage_subgraph.nodes[node_id]['blaster_file'])
                            new_blaster_file.path = f"{new_blaster_file.path}_{blaster_file_node_id_to_numerical_suffix_dict[node_id]}"
                            new_dock_file_lineage_subgraph.nodes[node_id]['blaster_file'] = new_blaster_file
                        dock_file_lineage_subgraph = new_dock_file_lineage_subgraph

                        #
                        for u, v, data in dock_file_lineage_subgraph.edges(data=True):
                            step_hash_to_edges_dict[data["step_hash"]].append((u, v))

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
                            step_hash_to_step_class_instance_dict[step_hash] = step_class(**filter_kwargs_for_callable(kwargs, step_class))

                            #
                            for parent_node, child_node in edges:
                                dock_file_lineage_subgraph.get_edge_data(parent_node, child_node)[
                                    "step_instance"
                                ] = step_hash_to_step_class_instance_dict[step_hash]

                        #
                        for u, v, data in dock_file_lineage_subgraph.edges(data=True):
                            if graph.has_edge(u, v):
                                for attr in ['parameter', 'blaster_file']:
                                    for n in [u, v]:
                                        if dock_file_lineage_subgraph.nodes[n].get(attr) is not None:
                                            if dock_file_lineage_subgraph.nodes[n].get(attr) != graph.nodes[n].get(attr):
                                                raise Exception(f"`dock_file_lineage_subgraph` and `graph` have nodes with ID `{n}` in common but possess unequal attribute `{attr}`: {dock_file_lineage_subgraph.nodes[n].get(attr)} vs. {graph.nodes[n].get(attr)}")
                            if graph.hash_node(v):
                                for pred in graph.predecessors(v):
                                    if graph.nodes[pred].get(u_node_type) is not None:
                                        if data[u_node_type] == graph.nodes[pred][u_node_type]:
                                            continue
                                        if data.get('parameter') is not None:
                                            if data['parameter'].name == graph.nodes[pred]['parameter'].name:
                                                raise Exception(f"Nodes with ID `{v}` in common in `dock_file_lineage_subgraph` and `graph` have different parent parameter nodes: {data[u_node_type]} vs. {graph.nodes[pred][u_node_type]}")
                                        elif data.get('blaster_file') is not None:
                                            if data['blaster_file'].identifier == graph.nodes[pred]['blaster_file'].identifier:
                                                raise Exception(f"Nodes with ID `{v}` in common in `dock_file_lineage_subgraph` and `graph` have different parent blaster file nodes: {data[u_node_type]} vs. {graph.nodes[pred][u_node_type]}")
                                        else:
                                            raise Exception(f"Unrecognized node type for `{u}`: {data}")

                        #
                        graph = nx.compose(graph, dock_file_lineage_subgraph)

                #
                partial_dock_file_nodes_combination_dicts.append(partial_dock_file_nodes_combination_dict)

        #
        dc_kwargs_so_far = []
        if last_component_docking_configurations:
            if partial_dock_file_nodes_combination_dicts:
                for last_component_dc, partial_dock_file_nodes_combination_dict in itertools.product(last_component_docking_configurations, partial_dock_file_nodes_combination_dicts):
                    dock_file_coordinates_kwargs = {  # complement + complement = complete
                        **{identifier: DockFileCoordinate(
                            component_id=self.component_id,
                            file_name=graph.nodes[node_id]['blaster_file'].name,
                            node_id=node_id,
                        ) for identifier, node_id in partial_dock_file_nodes_combination_dict.items()},
                        **{identifier: coord for identifier, coord in asdict(last_component_dc.dock_file_coordinates).items() if identifier not in partial_dock_file_nodes_combination_dict},
                    }
                    dock_file_coordinates = DockFileCoordinates(**dock_file_coordinates_kwargs)
                    partial_dc_kwargs = {
                        'dock_file_coordinates': dock_file_coordinates,
                        'dock_files_generation_flat_param_dict': self._get_dock_files_generation_flat_param_dict(graph, dock_file_coordinates),
                    }
                    dc_kwargs_so_far.append(partial_dc_kwargs)
            else:
                for last_component_dc in last_component_docking_configurations:
                    partial_dc_kwargs = {
                        'dock_file_coordinates': deepcopy(last_component_dc.dock_file_coordinates),
                        'dock_files_generation_flat_param_dict': self._get_dock_files_generation_flat_param_dict(graph, last_component_dc.dock_file_coordinates),
                    }
                    dc_kwargs_so_far.append(partial_dc_kwargs)
        else:
            for partial_dock_file_nodes_combination_dict in partial_dock_file_nodes_combination_dicts:
                dock_file_coordinates_kwargs = {
                    **{identifier: DockFileCoordinate(
                        component_id=self.component_id,
                        file_name=graph.nodes[node_id]['blaster_file'].name,
                        node_id=node_id,
                    ) for identifier, node_id in partial_dock_file_nodes_combination_dict.items()},
                }
                dock_file_coordinates = DockFileCoordinates(**dock_file_coordinates_kwargs)
                partial_dc_kwargs = {
                    'dock_file_coordinates': dock_file_coordinates,
                    'dock_files_generation_flat_param_dict': self._get_dock_files_generation_flat_param_dict(graph, dock_file_coordinates),
                }
                dc_kwargs_so_far.append(partial_dc_kwargs)
        logger.debug(f"Number of partial docking configurations after dock files generation specification: {len(dc_kwargs_so_far)}")

        #
        dc_kwargs_so_far = self._get_unique_partial_docking_configuration_kwargs_sorted(dc_kwargs_so_far)
        logger.debug(f"Number of unique partial docking configurations after dock files generation specification: {len(dc_kwargs_so_far)}")

        # matching spheres perturbation
        sorted_unique_matching_spheres_file_nodes = sorted(list(set([partial_dc_kwargs['dock_file_coordinates'].matching_spheres_file.node_id for partial_dc_kwargs in dc_kwargs_so_far])))
        new_dc_kwargs_so_far = []
        num_files_perturbed_so_far = 0
        for dock_files_modification_flat_param_dict in sorted_dock_files_modification_flat_param_dicts:
            if dock_files_modification_flat_param_dict[
                "matching_spheres_perturbation.use"
            ].value:
                #
                matching_spheres_node_to_perturbed_nodes_dict = collections.defaultdict(list)
                for i in range(
                    int(
                        dock_files_modification_flat_param_dict[
                            "matching_spheres_perturbation.num_samples_per_matching_spheres_file"
                        ].value
                    )
                ):
                    for matching_spheres_file_node in sorted_unique_matching_spheres_file_nodes:
                        #
                        matching_spheres_blaster_file = graph.nodes[matching_spheres_file_node]['blaster_file']
                        max_deviation_angstroms = float(
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
                            f"{BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT[matching_spheres_blaster_file.identifier]}_p{num_files_perturbed_so_far+1}"  # 'p' for perturbed
                        )
                        perturbed_matching_spheres_file = BlasterFile(perturbed_matching_spheres_file_path, identifier="matching_spheres_file")
                        step = MatchingSpheresPerturbationStep(
                            self.working_dir,
                            matching_spheres_infile=matching_spheres_blaster_file,
                            perturbed_matching_spheres_outfile=perturbed_matching_spheres_file,
                            max_deviation_angstroms_parameter=max_deviation_angstroms_parameter,
                        )
                        num_files_perturbed_so_far += 1

                        # get step hash from infile hashes, step dir, parameters, and outfiles
                        step_hash = DockoptStep._get_step_hash(self.component_id, step)

                        # add outfile node
                        outfile, = list(step.outfiles._asdict().values())
                        outfile_hash = DockoptStep._get_outfile_hash(self.component_id, outfile, step_hash)
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
                            matching_spheres_file_node,
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
                        matching_spheres_node_to_perturbed_nodes_dict[matching_spheres_file_node].append(outfile_hash)

                #
                for partial_dc_kwargs in dc_kwargs_so_far:
                    dock_file_coordinates = partial_dc_kwargs['dock_file_coordinates']
                    for perturbed_file_node_id in matching_spheres_node_to_perturbed_nodes_dict[dock_file_coordinates.matching_spheres_file.node_id]:
                        new_partial_dc_kwargs = deepcopy(partial_dc_kwargs)
                        new_partial_dc_kwargs['dock_file_coordinates'].matching_spheres_file = DockFileCoordinate(
                            component_id=self.component_id,
                            file_name=graph.nodes[perturbed_file_node_id]['blaster_file'].name,
                            node_id=perturbed_file_node_id,
                        )
                        new_partial_dc_kwargs['dock_files_modification_flat_param_dict'] = dock_files_modification_flat_param_dict
                        new_dc_kwargs_so_far.append(new_partial_dc_kwargs)
            else:
                temp_dc_kwargs_so_far = deepcopy(dc_kwargs_so_far)
                for partial_dc_kwargs in temp_dc_kwargs_so_far:
                    partial_dc_kwargs['dock_files_modification_flat_param_dict'] = dock_files_modification_flat_param_dict
                new_dc_kwargs_so_far += temp_dc_kwargs_so_far
        logger.debug(f"Number of partial docking configurations after dock files modification specification: {len(new_dc_kwargs_so_far)}")

        #
        dc_kwargs_so_far = self._get_unique_partial_docking_configuration_kwargs_sorted(new_dc_kwargs_so_far)
        logger.debug(f"Number of unique partial docking configurations after dock files modification specification: {len(dc_kwargs_so_far)}")

        #
        new_dc_kwargs_so_far = []
        for i, (partial_dc_kwargs, custom_dock_executable, indock_file_generation_flat_param_dict) in enumerate(itertools.product(dc_kwargs_so_far, custom_dock_executables, sorted_indock_file_generation_flat_param_dicts)):
            configuration_num = i + 1
            new_partial_dc_kwargs = deepcopy(partial_dc_kwargs)
            new_partial_dc_kwargs = {
                'component_id': self.component_id,
                'configuration_num': configuration_num,
                'custom_dock_executable': custom_dock_executable,
                'dock_files_generation_flat_param_dict': partial_dc_kwargs['dock_files_generation_flat_param_dict'],
                'dock_files_modification_flat_param_dict': partial_dc_kwargs['dock_files_modification_flat_param_dict'],
                'indock_file_generation_flat_param_dict': indock_file_generation_flat_param_dict,
                'dock_file_coordinates': partial_dc_kwargs['dock_file_coordinates'],
                'indock_file_coordinate': IndockFileCoordinate(
                    component_id=self.component_id,
                    file_name=f"{INDOCK_FILE_NAME}_{configuration_num}",
                ),
            }
            new_dc_kwargs_so_far.append(new_partial_dc_kwargs)
        logger.debug(f"Number of partial docking configurations after indock file generation specification: {len(new_partial_dc_kwargs)}")

        #
        all_dc_kwargs = self._get_unique_partial_docking_configuration_kwargs_sorted(new_dc_kwargs_so_far)
        logger.debug(f"Number of unique partial docking configurations after indock file generation specification: {len(all_dc_kwargs)}")

        #
        self.docking_configurations = [DockingConfiguration(**dc_kwargs) for dc_kwargs in all_dc_kwargs]
        logger.info(f"Number of unique docking configurations: {len(self.docking_configurations)}")

        # validate that there are no cycles (i.e. that it is a directed acyclic graph)
        if not nx.is_directed_acyclic_graph(graph):
            raise Exception("Cycle found in graph!")

        #
        self.graph = graph
        logger.debug(
            f"Graph initialized with:\n\tNodes: {self.graph.nodes}\n\tEdges: {self.graph.edges}"
        )

    def _get_unique_partial_docking_configuration_kwargs_sorted(self, dc_kwargs_list: List[dict]) -> List[dict]:
        logger.debug(f"Getting unique partial docking configurations (sorted). # before: {len(dc_kwargs_list)}")
        new_dc_kwargs = []
        hashes = []
        for dc_kwargs in dc_kwargs_list:
            hash = DockingConfiguration.get_hexdigest_of_persistent_md5_hash_of_docking_configuration_kwargs(dc_kwargs, partial_okay=True)
            if hash not in hashes:
                new_dc_kwargs.append(dc_kwargs)
                hashes.append(hash)

        #
        new_dc_kwargs_sorted, hashes_sorted = zip(*sorted(zip(new_dc_kwargs, hashes), key=lambda x: x[1]))
        logger.debug(f"# after: {len(new_dc_kwargs_sorted)}")

        return new_dc_kwargs_sorted

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.DataFrame:
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
        for dc in self.docking_configurations:
            # make dock files
            for dock_file_identifier in DOCK_FILE_IDENTIFIERS:
                self._run_unrun_steps_needed_to_create_this_blaster_file_node(
                    getattr(dc.dock_file_coordinates, dock_file_identifier).node_id, self.graph
                )

            # make indock file now that dock files exist
            indock_file = dc.get_indock_file(self.pipeline_dir.path)
            indock_file.write(dc.get_dock_files(self.pipeline_dir.path), dc.indock_file_generation_flat_param_dict)
        logger.info("done.")

        #
        step_id_file_path = os.path.join(self.component_dir.path, "step_id")
        if File.file_exists(step_id_file_path):
            with open(step_id_file_path, "r") as f:
                step_id, = tuple([line.strip() for line in f.readlines()])
                try:
                    _ = uuid.UUID(step_id)
                except ValueError:
                    raise Exception("step id loaded from step_id_file_path is not a valid UUID.")
        else:
            step_id = str(uuid.uuid4())
            with open(step_id_file_path, "w") as f:
                f.write(f"{step_id}\n")

        #
        array_job_docking_configurations_file_path = os.path.join(self.component_dir.path, "array_job_docking_configurations.txt")
        with open(array_job_docking_configurations_file_path, 'w') as f:
            for dc in self.docking_configurations:
                dock_files = dc.get_dock_files(self.pipeline_dir.path)
                dockfile_paths_str = " ".join([getattr(dock_files, field.name).path for field in fields(dock_files)])
                indock_file_path_str = dc.get_indock_file(self.pipeline_dir.path).path
                f.write(f"{dc.configuration_num} {indock_file_path_str} {dockfile_paths_str} {dc.dock_executable_path}\n")

        #
        def get_actives_outdock_file_path_for_configuration_num(configuration_num: int) -> str:
            return os.path.join(self.retrodock_jobs_dir.path, 'actives', str(configuration_num), 'OUTDOCK.0')

        def get_decoys_outdock_file_path_for_configuration_num(configuration_num: int) -> str:
            return os.path.join(self.retrodock_jobs_dir.path, 'decoys', str(configuration_num), 'OUTDOCK.0')

        # submit retrodock jobs (one for actives, one for decoys)
        array_jobs = []
        for sub_dir_name, should_export_mol2, input_molecules_tgz_file_path in [
            ('actives', True, component_run_func_arg_set.actives_tgz_file_path),
            ('decoys', False, component_run_func_arg_set.decoys_tgz_file_path),
        ]:
            array_job = ArrayDockingJob(
                name=f"dockopt_step_{step_id}_{sub_dir_name}",
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
        datetime_queue_was_last_checked = datetime.min
        while len(docking_configurations_processing_queue) > 0:
            #
            docking_configuration = docking_configurations_processing_queue.pop(0)
            if any([not array_job.task_is_complete(str(docking_configuration.configuration_num)) for array_job in array_jobs]):  # one or both OUTDOCK files do not exist yet
                time.sleep(
                    0.01
                )  # sleep for a bit

                #
                datetime_now = datetime.now()
                if datetime_now > (datetime_queue_was_last_checked + timedelta(seconds=MIN_SECONDS_BETWEEN_QUEUE_CHECKS)):
                    datetime_queue_was_last_checked = datetime_now
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

            # make data dict for this configuration num
            data_dict = docking_configuration.to_dict()

            # get ROC and calculate enrichment score of this job's docking set-up
            if isinstance(self.criterion, EnrichmentScore):
                logger.debug("Calculating ROC and enrichment score...")
                booleans = df["is_active"]
                data_dict[self.criterion.name] = self.criterion.calculate(booleans)
                logger.debug("done.")

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
    def _get_dock_file_lineage_subgraph(graph: nx.DiGraph, dock_file_node_id: str) -> nx.DiGraph:
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
    def _get_infile_hash(component_id: str, infile: BlasterFile) -> str:
        return get_hexdigest_of_persistent_md5_hash_of_tuple((component_id, infile.original_file_in_working_dir.name))

    @staticmethod
    def _get_outfile_hash(component_id: str, outfile: BlasterFile, step_hash: str) -> str:
        return get_hexdigest_of_persistent_md5_hash_of_tuple((component_id, outfile.original_file_in_working_dir.name, step_hash))

    @staticmethod
    def _get_step_hash(component_id: str, step: BlasterStep) -> str:
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
    def _get_graph_from_all_steps_in_order(component_id: str, steps: List[BlasterStep]) -> nx.DiGraph:
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
    def _get_blaster_file_nodes(g: nx.DiGraph) -> str:
        return [node_id for node_id, node_data in g.nodes.items() if g.nodes[node_id].get("blaster_file")]

    @staticmethod
    def _get_blaster_file_node_with_blaster_file_identifier(
        blaster_file_identifier: str,
        g: nx.DiGraph,
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
        g: nx.DiGraph,
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
    def _get_dock_files_generation_flat_param_dict(graph: nx.DiGraph, dock_file_coordinates: DockFileCoordinates) -> dict:
        dock_file_node_ids = sorted([getattr(dock_file_coordinates, field.name).node_id for field in fields(dock_file_coordinates)])
        node_ids = [node_id for dock_file_node_id in dock_file_node_ids for node_id in nx.ancestors(graph, dock_file_node_id)]
        node_ids = list(set(node_ids))
        d = {}
        for node_id in node_ids:
            if graph.nodes[node_id].get('parameter'):
                parameter = graph.nodes[node_id]['parameter']
                d[parameter.name] = parameter.value
        logger.debug(f"dock files generation parameters dict derived from dock file nodes:\n\tnodes: {dock_file_node_ids}\n\tdict: {d}")

        return d

    @staticmethod
    def _run_unrun_steps_needed_to_create_this_blaster_file_node(
        blaster_file_node: str,
        g: nx.DiGraph,
    ) -> None:
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
        pipeline_dir_path: str,
        component_id: str,
        criterion: str,
        top_n: int,
        components: Iterable[dict],
        blaster_files_to_copy_in: Iterable[BlasterFile],
        last_component_completed: Union[PipelineComponent, None] = None,
    ) -> None:
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
            component_id=component_id,
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
            path=os.path.join(self.component_dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.graph = nx.DiGraph()

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.DataFrame:
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
            kwargs = filter_kwargs_for_callable({
                **parameters_manager.parameters_dict,
                'component_id': f"{self.component_id}.{component_num}",
                'pipeline_dir_path': self.pipeline_dir.path,
                'blaster_files_to_copy_in': self.blaster_files_to_copy_in,  # TODO: is this necessary?
                'last_component_completed': last_component_completed_in_sequence,
            }, component_class)
            component = component_class(**kwargs)
            component.run(component_run_func_arg_set)

            #
            df_component = component.load_results_dataframe()
            df = pd.concat([df, df_component], ignore_index=True)

            #
            last_component_completed_in_sequence = component

            # TODO: make sure this is memory efficient
            self.graph = nx.compose(self.graph, last_component_completed_in_sequence.graph)

        return df


class DockoptStepSequence(PipelineComponentSequence):
    def __init__(
        self,
        pipeline_dir_path: str,
        component_id: str,
        criterion: str,
        top_n: int,
        components: Iterable[dict],
        num_repetitions: int,
        max_iterations_with_no_improvement: int,
        inter_iteration_criterion: str,
        inter_iteration_top_n: int,
        blaster_files_to_copy_in: Iterable[BlasterFile],
        last_component_completed: Union[PipelineComponent, None] = None,
    ) -> None:
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
            component_id=component_id,
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
            path=os.path.join(self.component_dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.graph = nx.DiGraph()

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.DataFrame:
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
                pipeline_dir_path=self.pipeline_dir.path,
                criterion=component_criterion,
                top_n=component_top_n,
                components=self.components,
                blaster_files_to_copy_in=self.blaster_files_to_copy_in,
                last_component_completed=last_component_completed_in_sequence,
            )
            component.run(component_run_func_arg_set)

            #
            df_component = component.load_results_dataframe()
            df = pd.concat([df, df_component], ignore_index=True)

            #
            last_component_completed_in_sequence = component

            # TODO: make sure this is memory efficient
            self.graph = nx.compose(self.graph, last_component_completed_in_sequence.graph)

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


class DockoptPipeline(Pipeline):
    def __init__(
            self,
            pipeline_dir_path: str,
            criterion: str,
            top_n: int,
            components: Iterable[dict],
            blaster_files_to_copy_in: Iterable[BlasterFile],
            last_component_completed: Union[PipelineComponent, None] = None,
    ) -> None:
        super().__init__(
            pipeline_dir_path=pipeline_dir_path,
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
            path=os.path.join(self.pipeline_dir.path, BEST_RETRODOCK_JOBS_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        self.graph = nx.DiGraph()

    def run(self, component_run_func_arg_set: DockoptPipelineComponentRunFuncArgSet) -> pd.DataFrame:
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
                raise Exception(
                    f"Dict must have one of 'step' or 'sequence' as keys. Witnessed: {component_identifier_dict}")

            #
            component_parameters_dict = deepcopy(component_identifier_dict[component_identifier])
            parameters_manager = DockoptComponentParametersManager(
                parameters_dict=component_parameters_dict,
                last_component_completed=last_component_completed_in_sequence,
            )

            #
            kwargs = filter_kwargs_for_callable({
                **parameters_manager.parameters_dict,
                'component_id': str(component_num),
                'pipeline_dir_path': self.pipeline_dir.path,
                'blaster_files_to_copy_in': self.blaster_files_to_copy_in,  # TODO: is this necessary?
                'last_component_completed': last_component_completed_in_sequence,
            }, component_class)
            component = component_class(**kwargs)
            component.run(component_run_func_arg_set)

            #
            df_component = component.load_results_dataframe()
            df = pd.concat([df, df_component], ignore_index=True)

            #
            last_component_completed_in_sequence = component

            # TODO: make sure this is memory efficient
            self.graph = nx.compose(self.graph, last_component_completed_in_sequence.graph)

        return df
