import itertools
import os
import shutil
from functools import wraps
from dataclasses import fields
from copy import copy
import logging
import collections
import time
import random

import networkx as nx
import numpy as np
import pandas as pd
from dirhash import dirhash

from pydock3.util import CleanExit, get_dataclass_as_dict
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
from pydock3.blastermaster.util import WorkingDir, BlasterFile, DockFiles, BlasterFileNames
from pydock3.dockopt.roc import ROC
from pydock3.jobs import RetrodockJob
from pydock3.job_schedulers import SlurmJobScheduler, SGEJobScheduler
from pydock3.dockopt.report import generate_dockopt_job_report
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH
from pydock3.blastermaster.defaults import __file__ as DEFAULTS_INIT_FILE_PATH
from pydock3.blastermaster.programs.thinspheres.sph_lib import read_sph, write_sph


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#
SCHEDULER_NAME_TO_CLASS_DICT = {
    "sge": SGEJobScheduler,
    "slurm": SlurmJobScheduler,
}


class AbstractTargetsDAG(object):

    def __init__(self, config, blaster_files, working_dir):
        #
        ad_hoc_job_param_dict = config.job_param_dicts[0]

        #
        steps = get_blaster_steps(
            blaster_files=blaster_files,
            param_dict=ad_hoc_job_param_dict,
            working_dir=working_dir,
        )

        #
        g = nx.DiGraph()
        for step in steps:

            # add infile nodes
            for infile in step.infiles:
                g.add_node(infile.original_file_in_working_dir.name, target=infile.original_file_in_working_dir)

            # add outfile nodes
            for outfile in step.outfiles:
                g.add_node(outfile.original_file_in_working_dir.name, target=outfile.original_file_in_working_dir)

            # add parameter nodes
            for parameter in step.parameters:
                g.add_node(parameter.name, parameter=parameter)

            #
            infiles_dict_items_list = step.infiles._asdict().items()
            outfiles_dict_items_list = step.outfiles._asdict().items()
            parameters_dict_items_list = step.parameters._asdict().items()

            # connect each infile node to every outfile node
            for (infile_step_var_name, infile), (outfile_step_var_name, outfile) in itertools.product(infiles_dict_items_list, outfiles_dict_items_list):
                g.add_edge(
                    infile.original_file_in_working_dir.name,
                    outfile.original_file_in_working_dir.name,
                    step=step,
                    parent_node_step_var_name=infile_step_var_name,
                    child_node_step_var_name=outfile_step_var_name,
                )

            # connect each parameter node to every outfile nodes
            for (parameter_step_var_name, parameter), (outfile_step_var_name, outfile) in itertools.product(parameters_dict_items_list, outfiles_dict_items_list):
                g.add_edge(
                    parameter.name,
                    outfile.original_file_in_working_dir.name,
                    step=step,
                    parent_node_step_var_name=parameter_step_var_name,
                    child_node_step_var_name=outfile_step_var_name,
                )

        # validate that there are no cycles (i.e. that it is a DAG)
        if not nx.is_directed_acyclic_graph(g):
            raise Exception("Cycle found in blaster targets DAG!")

        self.g = g

        logger.debug(f"AbstractTargetsDAG initialized with:\n\tNodes: {self.g.nodes}\n\tEdges: {self.g.edges}")

        #
        self.dock_files = blaster_files.dock_files


class FullTargetsDAG(object):
    def __init__(self, abstract_targets_dag, config, working_dir):
        #
        self.g = nx.DiGraph()
        config = config
        working_dir = working_dir

        #
        abstract_nodes_processed = []

        #
        def process_abstract_node(abstract_node, abstract_targets_dag):
            #
            if abstract_node in abstract_nodes_processed:
                return

            #
            if abstract_targets_dag.g.in_degree(abstract_node) > 0:
                predecessors = list(abstract_targets_dag.g.predecessors(abstract_node))
                unprocessed_predecessors = [pred for pred in predecessors if pred not in abstract_nodes_processed]
                for unprocessed_pred in unprocessed_predecessors:  # all predecessor nodes need to be processed first
                    process_abstract_node(unprocessed_pred, abstract_targets_dag)

            #
            abstract_node_data_dict = abstract_targets_dag.g.nodes[abstract_node]
            if "target" in abstract_node_data_dict:
                abstract_target = abstract_node_data_dict["target"]

                # yes predecessors, add node and preceding edges for every combination
                if abstract_targets_dag.g.in_degree(abstract_node) > 0:
                    abstract_predecessors = list(abstract_targets_dag.g.predecessors(abstract_node))  # all predecessors have been processed
                    abstract_pred_node_to_corresponding_full_dag_pred_nodes_dict = {}
                    for abstract_pred in abstract_predecessors:
                        abstract_pred_node_to_corresponding_full_dag_pred_nodes_dict[abstract_pred] = [full_dag_pred_node for full_dag_pred_node, data in self.g.nodes(data=True) if data.get("abstract_node_name") == abstract_pred]
                    items = list(abstract_pred_node_to_corresponding_full_dag_pred_nodes_dict.items())
                    abstract_pred_nodes = [x[0] for x in items]
                    full_dag_pred_nodes_lists = [x[1] for x in items]
                    for i, full_dag_pred_nodes in enumerate(itertools.product(*full_dag_pred_nodes_lists)):
                        # skip if any pred nodes have different ancestor nodes with same abstract_node_name
                        full_dag_ancestors = []
                        for full_dag_pred_node in full_dag_pred_nodes:
                            full_dag_ancestors.append(full_dag_pred_node)
                            full_dag_ancestors += list(nx.ancestors(self.g, full_dag_pred_node))
                        full_dag_ancestors_abstract_node_names = [self.g.nodes[full_dag_ancestor]['abstract_node_name'] for full_dag_ancestor in full_dag_ancestors]
                        if len(list(set(full_dag_ancestors))) != len(list(set(full_dag_ancestors_abstract_node_names))):
                            continue

                        #
                        full_dag_target_name = f"{abstract_target.name}_{i+1}"
                        full_dag_target = BlasterFile(
                            path=os.path.join(working_dir.path, full_dag_target_name),
                            src_file_path=abstract_target.src_file_path,
                        )
                        self.g.add_node(full_dag_target_name, target=full_dag_target, abstract_node_name=abstract_node)
                        for j, full_dag_pred_node in enumerate(full_dag_pred_nodes):
                            abstract_pred_node = abstract_pred_nodes[j]
                            abstract_edge = abstract_targets_dag.g.get_edge_data(abstract_pred_node, abstract_node)
                            self.g.add_edge(
                                full_dag_pred_node,
                                full_dag_target_name,
                                step_identity=f"{abstract_edge['step'].step_dir.name}_{i+1}",
                                step_type=abstract_edge['step'].step_dir.name,
                                step_class=abstract_edge['step'].__class__,
                                parent_node_step_var_name=abstract_edge['parent_node_step_var_name'],
                                child_node_step_var_name=abstract_edge['child_node_step_var_name'],
                            )

                # no predecessors, add single node
                else:
                    full_dag_target = BlasterFile(
                        path=os.path.join(working_dir.path, abstract_target.name),
                        src_file_path=abstract_target.src_file_path
                    )
                    self.g.add_node(abstract_node, target=full_dag_target, abstract_node_name=abstract_node)

            elif "parameter" in abstract_node_data_dict:
                #
                abstract_parameter = abstract_node_data_dict["parameter"]
                values = list(set([job_param_dict[abstract_parameter.name].value for job_param_dict in config.job_param_dicts]))
                parameters = [Parameter(
                    name=abstract_parameter.name,
                    value=value,
                ) for value in values]

                # add parameter nodes to full dag
                for i, parameter in enumerate(parameters):
                    self.g.add_node(f"{parameter.name}_{i+1}", parameter=parameter, abstract_node_name=abstract_node)

            else:
                raise Exception(f"Abstract targets DAG node found that is neither a parameter nor target: {abstract_node}")

            #
            abstract_nodes_processed.append(abstract_node)

        # get abstract start nodes
        abstract_start_nodes = []
        for node in abstract_targets_dag.g.nodes:
            if abstract_targets_dag.g.in_degree(node) == 0:
                abstract_start_nodes.append(node)

        # process start nodes
        for abstract_start_node in abstract_start_nodes:
            process_abstract_node(abstract_start_node, abstract_targets_dag)

        # process all other nodes
        while set(abstract_nodes_processed) != set(list(abstract_targets_dag.g.nodes)):
            for edge in nx.edge_bfs(abstract_targets_dag.g, abstract_start_nodes):
                abstract_in_node, abstract_out_node = edge
                process_abstract_node(abstract_out_node, abstract_targets_dag)

        # create step instance for each step identity with args taken from all relevant edges, attach step instance to each relevant edge
        step_identity_to_edges_dict = collections.defaultdict(list)
        step_identity_to_step_class_dict = {}
        for u, v, data in self.g.edges(data=True):
            step_identity_to_edges_dict[data["step_identity"]].append((u, v))
            step_identity_to_step_class_dict[data["step_identity"]] = data["step_class"]  # TODO: fix this lazy implementation

        for step_identity, edges in step_identity_to_edges_dict.items():
            kwargs = {"step_dir": Dir(path=os.path.join(working_dir.path, step_identity))}
            for (parent_node, child_node) in edges:
                edge_data_dict = self.g.get_edge_data(parent_node, child_node)
                parent_node_data_dict = self.g.nodes[parent_node]
                child_node_data_dict = self.g.nodes[child_node]
                parent_node_step_var_name = edge_data_dict["parent_node_step_var_name"]
                child_node_step_var_name = edge_data_dict["child_node_step_var_name"]
                if "target" in parent_node_data_dict:
                    kwargs[parent_node_step_var_name] = parent_node_data_dict["target"]
                if "parameter" in parent_node_data_dict:
                    kwargs[parent_node_step_var_name] = parent_node_data_dict["parameter"]
                if "target" in child_node_data_dict:
                    kwargs[child_node_step_var_name] = child_node_data_dict["target"]
                if "parameter" in child_node_data_dict:
                    kwargs[child_node_step_var_name] = child_node_data_dict["parameter"]

            step = step_identity_to_step_class_dict[step_identity](**kwargs)

            for parent_node, child_node in edges:
                self.g.get_edge_data(parent_node, child_node)["step"] = step

        def get_ancestor_parameters(nodes):
            #
            ancestor_parameter_nodes = []
            for node in nodes:
                for ancestor_node in nx.ancestors(self.g, node):
                    if "parameter" in self.g.nodes[ancestor_node]:
                        ancestor_parameter_nodes.append(ancestor_node)

            #
            ancestor_parameter_nodes = list(set(ancestor_parameter_nodes))

            #
            ancestor_parameters = [self.g.nodes[ancestor_parameter_node]["parameter"] for ancestor_parameter_node in ancestor_parameter_nodes]

            return ancestor_parameters

        #
        abstract_dag_dock_file_node_name_to_dock_file_field_name_dict = {getattr(abstract_targets_dag.dock_files, dock_file_field.name).name: dock_file_field.name for dock_file_field in fields(abstract_targets_dag.dock_files)}
        abstract_dag_dock_file_node_name_to_full_dag_dock_file_nodes_dict = {abstract_dag_dock_file_node_name: [full_dag_node for full_dag_node in self.g.nodes if self.g.nodes[full_dag_node].get("abstract_node_name") == abstract_dag_dock_file_node_name] for abstract_dag_dock_file_node_name in abstract_dag_dock_file_node_name_to_dock_file_field_name_dict}
        abstract_dag_dock_file_node_names, full_dag_dock_file_nodes_lists = zip(*abstract_dag_dock_file_node_name_to_full_dag_dock_file_nodes_dict.items())
        dock_files_combinations = []
        input_parameters_combinations = []
        for full_dag_dock_file_nodes in itertools.product(*full_dag_dock_file_nodes_lists):
            if not self.nodes_have_consistent_lineages(full_dag_dock_file_nodes):  # skip dockfile node combinations that have inconsistent lineages
                continue
            kwargs = {abstract_dag_dock_file_node_name_to_dock_file_field_name_dict[abstract_dag_dock_file_node_name] : self.g.nodes[full_dag_dock_file_node]['target'] for abstract_dag_dock_file_node_name, full_dag_dock_file_node in zip(abstract_dag_dock_file_node_names, full_dag_dock_file_nodes)}
            dock_files = DockFiles(**kwargs)
            dock_files_combinations.append(dock_files)

            #
            ancestor_parameters = get_ancestor_parameters(full_dag_dock_file_nodes)
            input_parameters_combinations.append(ancestor_parameters)

        #
        self.dock_files_combinations = dock_files_combinations
        self.input_parameters_combinations = input_parameters_combinations

    def get_start_nodes(self):
        start_nodes = []
        for node in self.g.nodes:
            if self.g.in_degree(node) == 0:
                start_nodes.append(node)
        return start_nodes

    def get_edges_in_lineage(self, node):
        edges = []
        for pred in self.g.predecessors(node):
            edges.append((pred, node))
        for ancestor in nx.ancestors(self.g, node):
            for ancestor_pred in self.g.predecessors(ancestor):
                edges.append((ancestor_pred, ancestor))
        return list(set(edges))

    def nodes_have_consistent_lineages(self, nodes):
        edges_in_lineages = list(set([edge for node in nodes for edge in self.get_edges_in_lineage(node)]))
        return len(list(set([self.g.edges[edge]['step_identity'] for edge in edges_in_lineages]))) == len(list(set([self.g.edges[edge]['step_type'] for edge in edges_in_lineages])))


class Dockopt(object):

    JOB_DIR_NAME = "dockopt_job"
    CONFIG_FILE_NAME = "dockopt_config.yaml"
    ACTIVES_TGZ_FILE_NAME = "actives.tgz"
    DECOYS_TGZ_FILE_NAME = "decoys.tgz"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.dirname(DOCKOPT_INIT_FILE_PATH), "default_dockopt_config.yaml")
    WORKING_DIR_NAME = "working"
    RETRODOCK_JOBS_DIR_NAME = "retrodock_jobs"
    DEFAULT_FILES_DIR_PATH = os.path.dirname(DEFAULTS_INIT_FILE_PATH)

    def __init__(self):
        pass

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
        job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create working dir & copy in blaster files
        blaster_file_names = list(get_dataclass_as_dict(BlasterFileNames()).values())
        backup_blaster_file_paths = [os.path.join(self.DEFAULT_FILES_DIR_PATH, blaster_file_name) for blaster_file_name in blaster_file_names]
        blaster_file_names_in_cwd = [f for f in blaster_file_names if os.path.isfile(f)]
        files_to_copy_str = '\n\t'.join(blaster_file_names_in_cwd)
        if blaster_file_names_in_cwd:
            logger.info(f"Copying the following files from current directory into job working directory:\n\t{files_to_copy_str}")
        else:
            logger.info(f"No blaster files detected in current working directory. Be sure to add them manually before running the job.")
        working_dir = WorkingDir(path=os.path.join(job_dir.path, self.WORKING_DIR_NAME), create=True, reset=False, files_to_copy_in=blaster_file_names_in_cwd, backup_files_to_copy_in=backup_blaster_file_paths)

        # copy in actives and decoys TGZ files
        tgz_files = [self.ACTIVES_TGZ_FILE_NAME, self.DECOYS_TGZ_FILE_NAME]
        tgz_file_names_in_cwd = [f for f in tgz_files if os.path.isfile(f)]
        tgz_file_names_not_in_cwd = [f for f in tgz_files if not os.path.isfile(f)]
        if tgz_file_names_in_cwd:
            files_to_copy_str = '\n\t'.join(tgz_file_names_in_cwd)
            logger.info(f"Copying the following files from current directory into job directory:\n\t{files_to_copy_str}")
            for tgz_file_name in tgz_file_names_in_cwd:
                job_dir.copy_in_file(tgz_file_name)
        if tgz_file_names_not_in_cwd:
            files_missing_str = '\n\t'.join(tgz_file_names_not_in_cwd)
            logger.info(f"The following required files were not found in current working directory. Be sure to add them manually to the job directory before running the job.\n\t{files_missing_str}")

        # create retrodock jobs dir
        retrodock_jobs_dir = Dir(path=os.path.join(job_dir.path, self.RETRODOCK_JOBS_DIR_NAME), create=True, reset=False)

        # write fresh config file from default file
        save_path = os.path.join(job_dir.path, self.CONFIG_FILE_NAME)
        DockoptParametersConfiguration.write_config_file(save_path, self.DEFAULT_CONFIG_FILE_PATH, overwrite=overwrite)

    @handle_run_func.__get__(0)
    def run(self,
            scheduler,
            job_dir_path=".",
            config_file_path=None,
            actives_tgz_file_path=None,
            decoys_tgz_file_path=None,
            retrodock_job_max_reattempts=0,
            retrodock_job_timeout_minutes=None,
            ):
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
            logger.error("Actives TGZ file and/or decoys TGZ file not found. Did you put them in the job directory?\nNote: if you do not have actives and decoys, please use blastermaster instead of dockopt.")
            return
        if scheduler not in SCHEDULER_NAME_TO_CLASS_DICT:
            logger.error(f"scheduler flag must be one of: {list(SCHEDULER_NAME_TO_CLASS_DICT.keys())}")
            return

        #
        try:
            scheduler = SCHEDULER_NAME_TO_CLASS_DICT[scheduler]()
        except KeyError:
            logger.error(f"The following environmental variables are required to use the {scheduler} job scheduler: {SCHEDULER_NAME_TO_CLASS_DICT[scheduler].REQUIRED_ENV_VAR_NAMES}")
            return

        #
        try:
            TEMP_STORAGE_PATH = os.environ["TEMP_STORAGE_PATH"]
        except KeyError:
            logger.error("The following environmental variables are required to submit retrodock jobs: TEMP_STORAGE_PATH")
            return

        #
        job_dir = Dir(job_dir_path, create=True, reset=False)
        working_dir = WorkingDir(path=os.path.join(job_dir.path, self.WORKING_DIR_NAME), create=True, reset=False)
        retrodock_jobs_dir = Dir(path=os.path.join(job_dir.path, self.RETRODOCK_JOBS_DIR_NAME), create=True, reset=False)

        #
        logger.info("Loading config file...")
        config = DockoptParametersConfiguration(config_file_path)
        logger.info("done.")

        #
        config_params_str = '\n'.join(
            [f"{param_name}: {param.value}" for param_name, param in config.param_dict.items()])
        logger.info(f"Parameters:\n{config_params_str}")

        #
        if actives_tgz_file_path is not None:
            actives_tgz_file = File(path=actives_tgz_file_path)
        else:
            actives_tgz_file = None
        if actives_tgz_file_path is not None:
            decoys_tgz_file = File(path=decoys_tgz_file_path)
        else:
            decoys_tgz_file = None

        #
        blaster_files = BlasterFiles(working_dir=working_dir)

        #
        logger.info("Loading blaster targets graph...")
        abstract_targets_dag = AbstractTargetsDAG(config=config, blaster_files=blaster_files, working_dir=working_dir)
        full_targets_dag = FullTargetsDAG(abstract_targets_dag=abstract_targets_dag, config=config, working_dir=working_dir)
        logger.info("done.")
        
        def run_edge_step(edge):
            if full_targets_dag.g.get_edge_data(*edge).get("step").is_done:
                return

            _, child_node = edge
            all_parent_nodes = full_targets_dag.g.predecessors(child_node)
            for parent_node in all_parent_nodes:
                parent_parent_nodes = full_targets_dag.g.predecessors(parent_node)
                for parent_parent_node in parent_parent_nodes:
                    run_edge_step((parent_parent_node, parent_node))

            logger.debug(f"Running step of full dag edge: {edge}")
            full_targets_dag.g.get_edge_data(*edge).get("step").run()

        # run edge steps
        logger.info("Running blaster steps / validating blaster files")
        for edge in nx.edge_bfs(full_targets_dag.g, full_targets_dag.get_start_nodes()):
            run_edge_step(edge)

        # matching spheres perturbation
        dock_files_combinations_for_retro_docking = None
        input_parameters_combinations_for_retro_docking = None
        if config.param_dict["matching_spheres_perturbation.use"].value:
            #
            unperturbed_file_name_to_perturbed_file_names_dict = collections.defaultdict(list)
            for node_name, node_data in full_targets_dag.g.nodes(data=True):
                if node_data.get("abstract_node_name") == blaster_files.matching_spheres_file.name:
                    spheres = read_sph(os.path.join(working_dir.path, node_name), chosen_cluster='A', color='A')

                    #
                    for i in range(int(config.param_dict["matching_spheres_perturbation.num_samples_per_matching_spheres_file"].value)):

                        # perturb all spheres in file
                        new_spheres = []
                        for sphere in spheres:
                            new_sphere = copy(sphere)
                            perturbation_xyz = tuple([random.gauss(0.0, float(config.param_dict["matching_spheres_perturbation.standard_deviation_angstroms"].value)) for _ in range(3)])
                            new_sphere.X += perturbation_xyz[0]
                            new_sphere.Y += perturbation_xyz[1]
                            new_sphere.Z += perturbation_xyz[2]
                            new_spheres.append(new_sphere)

                        # write perturbed spheres to new matching spheres file
                        perturbed_file_name = f"{node_name}_{i+1}"
                        write_sph(os.path.join(working_dir.path, perturbed_file_name), new_spheres)

                        #
                        unperturbed_file_name_to_perturbed_file_names_dict[node_name].append(perturbed_file_name)

            #
            dock_files_combinations_for_retro_docking = []
            input_parameters_combinations_for_retro_docking = []
            for i, dock_files_combination in enumerate(full_targets_dag.dock_files_combinations):
                for perturbed_file_name in unperturbed_file_name_to_perturbed_file_names_dict[dock_files_combination.matching_spheres_file.name]:
                    new_dock_files_combination = copy(dock_files_combination)
                    new_dock_files_combination.matching_spheres_file = BlasterFile(os.path.join(working_dir.path, perturbed_file_name))
                    dock_files_combinations_for_retro_docking.append(new_dock_files_combination)
                    input_parameters_combinations_for_retro_docking.append(full_targets_dag.input_parameters_combinations[i])
        else:
            dock_files_combinations_for_retro_docking = full_targets_dag.dock_files_combinations
            input_parameters_combinations_for_retro_docking = full_targets_dag.input_parameters_combinations

        #
        matching_spheres_perturbation_param_dict = {key: value for key, value in config.param_dict.items() if key.startswith("matching_spheres_perturbation.")}

        #
        # TODO: there is an assumption here that none of the indock config params are used in any of the blaster steps; if this assumption ever breaks, this method will need to be replaced.
        job_param_dicts_indock_subset = [{key: value for key, value in job_param_dict.items() if key.startswith("indock.")} for job_param_dict in config.job_param_dicts]
        job_param_dicts_indock_subset = [dict(s) for s in set(frozenset(job_param_dict.items()) for job_param_dict in job_param_dicts_indock_subset)]  # get unique dicts

        #
        docking_configuration_info_combinations = list(itertools.product(zip(dock_files_combinations_for_retro_docking, input_parameters_combinations_for_retro_docking), job_param_dicts_indock_subset))

        # make indock file for each combination of (1) set of dock files and (2) job_param_dict_indock_subset
        logger.info("Making INDOCK files...")
        parameter_dicts = []
        docking_configurations = []
        for i, ((dock_files, input_parameters), job_param_dict_indock_subset) in enumerate(docking_configuration_info_combinations):
            # get full parameter dict
            parameter_dict = {p.name: p.value for p in input_parameters}
            parameter_dict.update(matching_spheres_perturbation_param_dict) # add matching spheres perturbation params
            parameter_dict.update(job_param_dict_indock_subset)  # add indock params

            # make indock file for each combination of dock files
            indock_file_name = f"{INDOCK_FILE_NAME}_{i+1}"
            indock_file = IndockFile(path=os.path.join(working_dir.path, indock_file_name))
            indock_file.write(dock_files, parameter_dict)

            #
            parameter_dicts.append(parameter_dict)
            docking_configurations.append((dock_files, indock_file))

        #
        all_docking_configuration_file_names = []
        for dock_files, indock_file in docking_configurations:
            dock_file_names = [getattr(dock_files, dock_file_field.name).name for dock_file_field in fields(dock_files)]
            all_docking_configuration_file_names += dock_file_names
            all_docking_configuration_file_names.append(indock_file.name)
        all_docking_configuration_file_names = list(set(all_docking_configuration_file_names))

        #
        dockopt_job_hash = dirhash(
            working_dir.path,
            "md5",
            match=all_docking_configuration_file_names,
        )
        with open(os.path.join(job_dir.path, "dockopt_job_hash.md5"), 'w') as f:
            f.write(f"{dockopt_job_hash}\n")

        # write actives tgz and decoys tgz file paths to actives_and_decoys.sdi
        logger.info("Writing actives_and_decoys.sdi file...")
        retrodock_input_sdi_file = File(path=os.path.join(job_dir.path, "actives_and_decoys.sdi"))
        with open(retrodock_input_sdi_file.path, 'w') as f:
            f.write(f"{actives_tgz_file.path}\n")
            f.write(f"{decoys_tgz_file.path}\n")
        logger.info("done")

        #
        retrodock_jobs = []
        retrodock_job_dirs = []
        retrodock_job_num_to_docking_configuration_file_names_dict = {}
        for i, (dock_files, indock_file) in enumerate(docking_configurations):
            #
            retro_dock_job_num = str(i+1)
            docking_configuration_file_names = [getattr(dock_files, dock_file_field.name).name for dock_file_field in fields(dock_files)] + [indock_file.name]
            retrodock_job_num_to_docking_configuration_file_names_dict[retro_dock_job_num] = docking_configuration_file_names

            #
            retrodock_job_dir = Dir(path=os.path.join(retrodock_jobs_dir.path, retro_dock_job_num), create=True)
            retrodock_job_output_dir = Dir(path=os.path.join(retrodock_job_dir.path, f"output"), create=True)

            #
            retrodock_job = RetrodockJob(
                name=f"dockopt_job_{dockopt_job_hash}_{retrodock_job_dir.name}",
                input_sdi_file=retrodock_input_sdi_file,
                dock_files=dock_files,
                indock_file=indock_file,
                output_dir=retrodock_job_output_dir,
                job_scheduler=scheduler,
                temp_storage_path=TEMP_STORAGE_PATH,
                max_reattempts=retrodock_job_max_reattempts,
            )
            retrodock_jobs.append(retrodock_job)
            retrodock_job_dirs.append(retrodock_job_dir)
        logger.debug("done")

        #
        def submit_retrodock_job(retrodock_job, skip_if_complete):
            logger.info(f"Submitting docking job for {retrodock_job.name}...")
            proc = retrodock_job.run(job_timeout_minutes=retrodock_job_timeout_minutes, skip_if_complete=skip_if_complete)
            if proc is None:
                logger.info(
                    f"Skipping docking job submission for {retrodock_job.name} since all its OUTDOCK files already exist.\n")
            else:
                logger.debug(
                    f"Retrodock job submission system call returned: {proc}\n\nstdout:{proc.stdout}\n\nstderr:{proc.stderr}\n")
                if proc.stderr:
                    logger.info(f"Job submission failed due to error: {proc.stderr}\n")
                else:
                    logger.info("done.\n")

        # submit docking jobs
        for retrodock_job in retrodock_jobs:
            submit_retrodock_job(retrodock_job, skip_if_complete=True)
        
        # make a queue of tuples containing job-relevant data for processing
        RetrodockJobInfoTuple = collections.namedtuple("RetrodockJobInfoTuple", "job job_dir parameter_dict")
        retrodock_jobs_processing_queue = [RetrodockJobInfoTuple(retrodock_jobs[i], retrodock_job_dirs[i], parameter_dicts[i]) for i in range(len(retrodock_jobs))]

        # process results of docking jobs
        logger.info(f"Awaiting / processing retrodock job results ({len(retrodock_job_dirs)} jobs in total)")
        data_dicts = []
        while len(retrodock_jobs_processing_queue) > 0:
            #
            retrodock_job_info_tuple = retrodock_jobs_processing_queue.pop(0)
            retrodock_job, retrodock_job_dir, parameter_dict = retrodock_job_info_tuple

            #
            if retrodock_job.is_running:
                retrodock_jobs_processing_queue.append(retrodock_job_info_tuple)  # move job to back of queue
                time.sleep(1)  # sleep a bit while waiting for outdock file in order to avoid wasteful queue-cycling
                continue  # move on to next job in queue while job continues to run
            else:
                if not retrodock_job.is_complete:  # not all expected OUTDOCK files exist yet
                    # job must have timed out / failed
                    logger.warning(f"Job failure / time out witnessed for job: {retrodock_job.name}")
                    if retrodock_job.num_attempts > retrodock_job_max_reattempts:
                        logger.warning(f"Max job reattempts exhausted for job: {retrodock_job.name}")
                        continue  # move on to next job in queue without re-attempting failed job
                    submit_retrodock_job(retrodock_job, skip_if_complete=False)  # re-attempt job
                    retrodock_jobs_processing_queue.append(retrodock_job_info_tuple)  # move job to back of queue
                    continue  # move on to next job in queue while docking job runs

            #
            actives_outdock_file_path = os.path.join(retrodock_job.output_dir.path, "1", "OUTDOCK.0")
            decoys_outdock_file_path = os.path.join(retrodock_job.output_dir.path, "2", "OUTDOCK.0")

            # load outdock file and get dataframe
            actives_outdock_file = OutdockFile(actives_outdock_file_path)
            decoys_outdock_file = OutdockFile(decoys_outdock_file_path)
            try:
                actives_outdock_df = actives_outdock_file.get_dataframe()
                decoys_outdock_df = decoys_outdock_file.get_dataframe()
            except Exception as e:  # if outdock file failed to be parsed then re-attempt job
                logger.warning(f"Failed to parse outdock file(s) due to error: {e}")
                if retrodock_job.num_attempts > retrodock_job_max_reattempts:
                    logger.warning(f"Max job reattempts exhausted for job: {retrodock_job.name}")
                    continue  # move on to next job in queue without re-attempting failed job
                submit_retrodock_job(retrodock_job, skip_if_complete=False)  # re-attempt job
                retrodock_jobs_processing_queue.append(retrodock_job_info_tuple)  # move job to back of queue
                continue  # move on to next job in queue while docking job runs

            #
            logger.info(f"Docking job '{retrodock_job.name}' completed. Successfully loaded OUTDOCK file(s).")

            # set is_active column based on outdock file
            actives_outdock_df["is_active"] = [1 for _ in range(len(actives_outdock_df))]
            decoys_outdock_df["is_active"] = [0 for _ in range(len(decoys_outdock_df))]

            # build dataframe of docking results from outdock files
            df = pd.DataFrame()
            df = pd.concat([df, actives_outdock_df], ignore_index=True)
            df = pd.concat([df, decoys_outdock_df], ignore_index=True)

            # sort dataframe by total energy score
            df["Total"] = df["Total"].astype(float)
            df = df.sort_values(by=["Total", "is_active"], na_position="last", ignore_index=True)  # sorting secondarily by 'is_active' (0 or 1) ensures that decoys are ranked before actives in case they have the same exact score (pessimistic approach)
            df = df.drop_duplicates(subset=["db2_file_path"], keep="first", ignore_index=True)

            # make data dict for this job (will be used to make dataframe for results of all jobs)
            data_dict = copy(parameter_dict)
            data_dict['retrodock_job_num'] = retrodock_job_dir.name

            # get ROC and calculate enrichment score of this job's docking set-up
            logger.debug("Calculating ROC and enrichment score...")
            booleans = df["is_active"]
            indices = df['Total'].fillna(np.inf)  # unscored molecules are assumed to have worst possible score (pessimistic approach)
            roc = ROC(booleans, indices)
            data_dict["enrichment_score"] = roc.enrichment_score
            logger.debug("done.")
            
            # write ROC plot image
            roc_plot_image_path = os.path.join(retrodock_job_dir.path, "roc.png")
            roc.plot(save_path=roc_plot_image_path)

            # save data_dict for this job
            data_dicts.append(data_dict)
           
        # write jobs completion status
        num_jobs_completed = len([1 for retrodock_job in retrodock_jobs if retrodock_job.is_complete])
        logger.info(f"Finished {num_jobs_completed} out of {len(retrodock_jobs)} retrodock jobs.")
        if num_jobs_completed == 0:
            logger.error("All retrodock jobs failed. Something is wrong. Please check logs.")
            return

        # make dataframe of optimization job results
        df = pd.DataFrame(data=data_dicts)
        df = df.sort_values(by="enrichment_score", ascending=False, ignore_index=True)

        # save optimization job results dataframe to csv
        optimization_results_csv_file_path = os.path.join(job_dir.path, "dockopt_job_results.csv")
        logger.debug(f"Saving optimization job results to {optimization_results_csv_file_path}")
        df.to_csv(optimization_results_csv_file_path)

        # copy best job to output dir
        best_retrodock_job_dir = Dir(os.path.join(job_dir.path, "best_retrodock_job"))
        logger.debug(f"Copying dockfiles of best job results to {best_retrodock_job_dir.path}")
        if os.path.isdir(best_retrodock_job_dir.path):
            shutil.rmtree(best_retrodock_job_dir.path, ignore_errors=True)
        best_retrodock_job_num = df['retrodock_job_num'].iloc[0]
        shutil.copytree(os.path.join(retrodock_jobs_dir.path, best_retrodock_job_num), best_retrodock_job_dir.path)

        # copy docking configuration files to best job dir
        best_retrodock_job_dockfiles_dir = Dir(os.path.join(best_retrodock_job_dir.path, "dockfiles"), create=True)
        for file_name in retrodock_job_num_to_docking_configuration_file_names_dict[best_retrodock_job_num]:
            best_retrodock_job_dockfiles_dir.copy_in_file(os.path.join(working_dir.path, file_name))

        # generate report
        generate_dockopt_job_report(
            dockopt_job_dir_path=job_dir.path,
            pdf_path=os.path.join(job_dir.path, "dockopt_job_report.pdf"),
            opt_results_csv_file_name=optimization_results_csv_file_path,
            enrichment_metric="enrichment_score",
        )
