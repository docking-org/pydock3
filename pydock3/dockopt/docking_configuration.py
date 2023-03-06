from dataclasses import dataclass
import collections
import itertools

from pydock3.files import IndockFile
from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple
from pydock3.blastermaster.util import DOCK_FILE_IDENTIFIERS, DockFiles

#
DockFileNodesTuple = collections.namedtuple("DockFileNodesTuple", " ".join(DOCK_FILE_IDENTIFIERS))


@dataclass
class DockingConfiguration:
    component_id: str
    configuration_num: str
    dock_executable_path: str
    dock_files_generation_flat_param_dict: dict
    dock_files_modification_flat_param_dict: dict
    indock_file_generation_flat_param_dict: dict
    dock_file_nodes_tuple: DockFileNodesTuple
    dock_files: DockFiles
    indock_file: IndockFile

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
        return get_hexdigest_of_persistent_md5_hash_of_tuple(tuple([getattr(self.dock_file_nodes_tuple, dock_file_identifier) for dock_file_identifier in DOCK_FILE_IDENTIFIERS] + [parameters_dict_items_interleaved_sorted_by_key_tuple]))
