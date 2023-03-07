from dataclasses import dataclass
import collections
import itertools
import os

from pydock3.files import IndockFile
from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple
from pydock3.blastermaster.util import DOCK_FILE_IDENTIFIERS, DockFiles
from pydock3.dockopt.util import WORKING_DIR_NAME

#
DockFileNodesTuple = collections.namedtuple("DockFileNodesTuple", " ".join(DOCK_FILE_IDENTIFIERS))

#
DockFileCoordinate = collections.namedtuple("DockFileCoordinate", "component_id file_name node_id")
IndockFileCoordinate = collections.namedtuple("IndockFileCoordinate", "component_id file_name")

#
DockFileCoordinates = collections.namedtuple("DockFileCoordinates", " ".join(DOCK_FILE_IDENTIFIERS))


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

    @property
    def full_flat_parameters_dict(self):
        flat_param_dict = {}
        flat_param_dict["dock_executable_path"] = self.dock_executable_path
        flat_param_dict.update(
            {
                f"dock_files_generation.{key}": value
                for key, value in self.dock_files_generation_flat_param_dict.items()
            }
        )
        flat_param_dict.update(
            {
                f"dock_files_modification.{key}": value
                for key, value in self.dock_files_modification_flat_param_dict.items()
            }
        )
        flat_param_dict.update(
            {
                f"indock_file_generation.{key}": value
                for key, value in self.indock_file_generation_flat_param_dict.items()
            }
        )

        return flat_param_dict

    @property
    def hexdigest_of_persistent_md5_hash(self):
        parameters_dict_items_interleaved_sorted_by_key_tuple = tuple(
            itertools.chain.from_iterable(
                sorted(list(zip(*list(zip(*self.full_flat_parameters_dict.items())))), key=lambda x: x[0])
            )
        )

        # order matters, so use DOCK_FILE_IDENTIFIERS
        return get_hexdigest_of_persistent_md5_hash_of_tuple(tuple([coordinate.node_id for dock_file_identifier, coordinate in self.dock_file_coordinates._asdict().items()] + [parameters_dict_items_interleaved_sorted_by_key_tuple]))

    @property
    def to_dict(self):
        d = {
            'component_id': self.component_id,
            'configuration_num': self.configuration_num,
            **self.full_flat_parameters_dict,
        }
        for dock_file_identifier, dock_file_coordinate in self.dock_file_coordinates:
            identifier_str = f"dock_files.{dock_file_identifier}"
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
                field: d[f"dock_files.{dock_file_identifier}.{field}"]
                for field in DockFileCoordinate._fields
            }) for dock_file_identifier in DOCK_FILE_IDENTIFIERS
        })
        indock_file_coordinate = IndockFileCoordinate(**{field: d[f"indock_file.{field}"] for field in DockFileCoordinate._fields})
        dock_files_generation_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.dock_files_generation')}
        dock_files_modification_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.dock_files_modification')}
        indock_file_generation_flat_param_dict = {key: value for key, value in d.items() if key.startswith('parameters.indock_file_generation')}
        return DockingConfiguration(
            component_id=d['component_id'],
            configuration_num=d['configuration_num'],
            dock_executable_path=d['parameters.dock_executable_path'],
            dock_files_generation_flat_param_dict=dock_files_generation_flat_param_dict,
            dock_files_modification_flat_param_dict=dock_files_modification_flat_param_dict,
            indock_file_generation_flat_param_dict=indock_file_generation_flat_param_dict,
            dock_file_coordinates=dock_file_coordinates,
            indock_file_coordinate=indock_file_coordinate,
        )

    def get_dock_files(self):
        kwargs = {dock_file_identifier: os.path.join(coordinate.component_id, WORKING_DIR_NAME, coordinate.file_name) for dock_file_identifier, coordinate in self.dock_file_coordinates._asdict().items()}
        return DockFiles(**kwargs)

    def get_indock_file(self):
        return IndockFile(os.path.join(self.indock_file_coordinate.component_id, WORKING_DIR_NAME, self.indock_file_coordinate.file_name))
