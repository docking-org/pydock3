from dataclasses import dataclass, make_dataclass, fields, asdict
import itertools
import os
import hashlib

from pydock3.files import IndockFile
from pydock3.blastermaster.util import BlasterFile
from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple, filter_kwargs_for_callable
from pydock3.blastermaster.util import DOCK_FILE_IDENTIFIERS, DockFiles
from pydock3.dockopt.util import WORKING_DIR_NAME


#
DockFileCoordinate = make_dataclass("DockFileCoordinate", [
    ("component_id", str),
    ("file_name", str),
    ("node_id", str),
])
IndockFileCoordinate = make_dataclass("IndockFileCoordinate", [
    ("component_id", str),
    ("file_name", str),
])

#
DockFileCoordinates = make_dataclass("DockFileCoordinates", [(identifier, DockFileCoordinate) for identifier in DOCK_FILE_IDENTIFIERS])


@dataclass
class DockingConfiguration:
    component_id: str
    configuration_num: str
    dock_executable_path: str
    dock_files_generation_flat_param_dict: dict
    dock_files_modification_flat_param_dict: dict
    indock_file_generation_flat_param_dict: dict
    dock_file_coordinates: DockFileCoordinates
    indock_file_coordinate: IndockFileCoordinate

    @staticmethod
    def get_full_flat_parameters_dict(dock_executable_path, dock_files_generation_flat_param_dict, dock_files_modification_flat_param_dict, indock_file_generation_flat_param_dict):
        flat_param_dict = {}
        flat_param_dict["dock_executable_path"] = dock_executable_path
        flat_param_dict.update(
            {
                f"dock_files_generation.{key}": value
                for key, value in dock_files_generation_flat_param_dict.items()
            }
        )
        flat_param_dict.update(
            {
                f"dock_files_modification.{key}": value
                for key, value in dock_files_modification_flat_param_dict.items()
            }
        )
        flat_param_dict.update(
            {
                f"indock_file_generation.{key}": value
                for key, value in indock_file_generation_flat_param_dict.items()
            }
        )

        return flat_param_dict

    @property
    def full_flat_parameters_dict(self):
        return self.get_full_flat_parameters_dict(self.dock_executable_path, self.dock_files_generation_flat_param_dict, self.dock_files_modification_flat_param_dict, self.indock_file_generation_flat_param_dict)

    @staticmethod
    def get_hexdigest_of_persistent_md5_hash_of_docking_configuration_kwargs(dc_kwargs):
        #
        filtered_kwargs = filter_kwargs_for_callable(dc_kwargs, DockingConfiguration.get_full_flat_parameters_dict)
        full_flat_parameters_dict = DockingConfiguration.get_full_flat_parameters_dict(**filtered_kwargs)

        #
        flat_parameters_dict_no_dock_exec_path = {key: value for key, value in full_flat_parameters_dict.items() if key != "dock_executable_path"}

        #
        parameters_dict_items_interleaved_sorted_by_key_tuple = tuple(
            itertools.chain.from_iterable(
                sorted(list(zip(*list(zip(*flat_parameters_dict_no_dock_exec_path.items())))), key=lambda x: x[0])
            )
        )

        #
        dock_exec_hash_tuple = tuple(hashlib.md5(open(full_flat_parameters_dict["dock_executable_path"], 'rb').read()).hexdigest())

        #
        dock_file_nodes_tuple = tuple([getattr(dc_kwargs['dock_file_coordinates'], field.name).node_id for field in fields(dc_kwargs['dock_file_coordinates'])])

        return get_hexdigest_of_persistent_md5_hash_of_tuple(parameters_dict_items_interleaved_sorted_by_key_tuple + dock_exec_hash_tuple + dock_file_nodes_tuple)

    @property
    def hexdigest_of_persistent_md5_hash(self):
        return self.get_hexdigest_of_persistent_md5_hash_of_docking_configuration_kwargs({field.name: getattr(self, field.name) for field in fields(self)})

    def to_dict(self):
        d = {
            'component_id': self.component_id,
            'configuration_num': self.configuration_num,
            **{f"parameters.{key}": value for key, value in self.full_flat_parameters_dict.items()},
        }
        for field in fields(self.dock_file_coordinates):
            dock_file_coordinate = getattr(self.dock_file_coordinates, field.name)
            identifier_str = f"dock_files.{field.name}"
            d[f"{identifier_str}.component_id"] = dock_file_coordinate.component_id
            d[f"{identifier_str}.file_name"] = dock_file_coordinate.file_name
            d[f"{identifier_str}.node_id"] = dock_file_coordinate.node_id
        d[f"indock_file.component_id"] = self.indock_file_coordinate.component_id
        d[f"indock_file.file_name"] = self.indock_file_coordinate.file_name

        return d

    @staticmethod
    def from_dict(d):
        dock_file_coordinates = DockFileCoordinates(**{
            dock_file_identifier: DockFileCoordinate(**{
                field.name: str(d[f"dock_files.{dock_file_identifier}.{field.name}"])
                for field in fields(DockFileCoordinate)
            }) for dock_file_identifier in DOCK_FILE_IDENTIFIERS
        })
        indock_file_coordinate = IndockFileCoordinate(**{field.name: str(d[f"indock_file.{field.name}"]) for field in fields(IndockFileCoordinate)})
        dock_files_generation_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.dock_files_generation')}
        dock_files_modification_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.dock_files_modification')}
        indock_file_generation_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.indock_file_generation')}
        return DockingConfiguration(
            component_id=str(d['component_id']),
            configuration_num=d['configuration_num'],
            dock_executable_path=d['parameters.dock_executable_path'],
            dock_files_generation_flat_param_dict=dock_files_generation_flat_param_dict,
            dock_files_modification_flat_param_dict=dock_files_modification_flat_param_dict,
            indock_file_generation_flat_param_dict=indock_file_generation_flat_param_dict,
            dock_file_coordinates=dock_file_coordinates,
            indock_file_coordinate=indock_file_coordinate,
        )

    def get_dock_files(self, pipeline_dir_path):
        kwargs = {field.name: BlasterFile(os.path.join(pipeline_dir_path, *getattr(self.dock_file_coordinates, field.name).component_id.split('.'), WORKING_DIR_NAME, getattr(self.dock_file_coordinates, field.name).file_name), identifier=field.name) for field in fields(self.dock_file_coordinates)}
        return DockFiles(**kwargs)

    def get_indock_file(self, pipeline_dir_path):
        return IndockFile(os.path.join(pipeline_dir_path, *self.indock_file_coordinate.component_id.split('.'), WORKING_DIR_NAME, self.indock_file_coordinate.file_name))
