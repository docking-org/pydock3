import copy
from typing import List, Tuple, Union, Optional, Dict, Any, TextIO, Generator
from enum import Enum
import collections
import logging
import os
import shutil
import pathlib
from datetime import datetime
import tarfile
import gzip
import re
import uuid
import time

import numpy as np
import pandas as pd
from rdkit import Chem

from pydock3.util import validate_variable_type, system_call


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
INDOCK_FILE_NAME = "INDOCK"


def create_relative_symlink(target: str, link_name: str, target_is_directory: bool) -> None:
    # Convert target and link_name to absolute paths
    target_absolute = pathlib.Path(target).resolve()
    link_name_absolute = pathlib.Path(link_name).resolve()

    # Calculate the relative path from link_name to target
    target_relative = os.path.relpath(target_absolute, start=link_name_absolute.parent)

    # Create the symbolic link
    os.symlink(target_relative, link_name_absolute, target_is_directory=target_is_directory)


class FileSystemEntity(object):
    """E.g., file, directory, symlink"""

    def __init__(self, path: str, validate_existence: bool = False):
        self.path = path

        if validate_existence:
            if not os.path.exists(self.path):
                raise FileNotFoundError(f"Path `{self.path}` does not exist.")

    def __str__(self):
        return self.path

    def __repr__(self):
        return self.path

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        self.validate_path(path)
        self._path = os.path.abspath(path)

    @property
    def name(self):
        raise NotImplementedError

    @property
    def exists(self):
        raise NotImplementedError

    @property
    def validate_existence(self):
        raise NotImplementedError

    @staticmethod
    def validate_path(file_path):
        validate_variable_type(file_path, allowed_instance_types=(str,))
        return file_path  # TODO: think more about this


class Dir(FileSystemEntity):
    """#TODO"""

    def __init__(self, path: str, validate_existence: bool = False, create: bool = False, reset: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

        #
        if create:
            self.create(reset=reset)

    @property
    def name(self):
        return Dir.extract_dir_name_from_dir_path(self.path)

    @staticmethod
    def extract_dir_name_from_dir_path(dir_path: str):
        return os.path.basename(dir_path)

    @property
    def exists(self):
        return Dir.dir_exists(self.path)

    @property
    def validate_existence(self):
        if not self.exists:
            raise Exception(f"Dir {self.path} does not exist.")

    @staticmethod
    def dir_exists(dir_path: str):
        return os.path.isdir(dir_path)

    def create(self, reset: bool = False):
        """#TODO"""
        if reset:
            self.delete()
            pathlib.Path(self.path).mkdir(parents=True, exist_ok=True)
            logger.debug(f"Reset directory {self}.")
        else:
            if os.path.exists(self.path):
                logger.debug(
                    f"Tried to create directory {self} with reset=False but directory already exists."
                )
            else:
                pathlib.Path(self.path).mkdir(parents=True, exist_ok=True)
                logger.debug(f"Created directory {self}")

    def delete(self):
        self.delete_dir(self.path)

    def reset(self):
        self.create(reset=True)

    def copy_in_file(self, src_file_path: str, dst_file_name: Optional[str] = None, overwrite: bool = True):
        """#TODO"""
        File.validate_file_exists(src_file_path)

        if dst_file_name is None:
            dst_file_name = File.get_file_name_of_file(src_file_path)
        dst_file = File(path=os.path.join(self.path, dst_file_name))
        dst_file.copy_from(src_file_path=src_file_path, overwrite=overwrite)

        return dst_file

    @staticmethod
    def delete_dir(dir_path: str) -> None:
        if os.path.exists(dir_path):
            '''
            shutil.rmtree(dir_path, ignore_errors=True)
            while os.path.isdir(dir_path):  # TODO: hmmm. is this valid?
                pass
            '''
            system_call(f"rm -rf {dir_path}")
            logger.debug(f"Deleted directory `{dir_path}`.")

    @staticmethod
    def reset_directory_cache(dir_path: str) -> None:
        """Reset the cache of files in the directory so that we can check for new files without dealing with distributed file system issues."""
        os.scandir(dir_path)
        time.sleep(0.01)

    @staticmethod
    def validate_obj_is_dir(obj: Any) -> None:
        validate_variable_type(obj, allowed_instance_types=(Dir,))


class File(FileSystemEntity):
    """#TODO"""

    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    @property
    def name(self):
        return File.get_file_name_of_file(self.path)

    @property
    def datetime_last_modified(self):
        return self.get_datetime_file_was_last_modified(self.path)

    @staticmethod
    def get_file_name_of_file(file_path: str):
        File.validate_path(file_path)
        return os.path.basename(os.path.abspath(file_path))

    @staticmethod
    def get_dir_path_of_file(file_path: str):
        File.validate_path(file_path)
        return os.path.dirname(os.path.abspath(file_path))

    @property
    def exists(self):
        return File.file_exists(self.path)

    @property
    def validate_existence(self):
        if not self.exists:
            raise Exception(f"File {self.path} does not exist.")

    @property
    def is_empty(self):
        return self.file_is_empty(self.path)

    def copy_from(self, src_file_path: str, overwrite: bool = True):
        self.copy_file(
            src_file_path=src_file_path, dst_file_path=self.path, overwrite=overwrite
        )

    def delete(self):
        self.delete_file(self.path)

    def validate_is_not_empty(self):
        if self.is_empty:
            raise Exception(f"File is empty: {self}")

    def read_lines(self):
        return self.read_file_lines(self.path)

    @property
    def is_gzipped(self):
        return self.file_is_gzipped(self.path)

    @staticmethod
    def get_datetime_file_was_last_modified(file_path: str):
        File.validate_file_exists(file_path)

        datetime_last_modified = datetime.fromtimestamp(os.stat(file_path).st_mtime)
        logger.debug(f"File {file_path} was last modified at: {datetime_last_modified}")
        return datetime_last_modified

    @staticmethod
    def get_file_size(file_path: str):
        File.validate_file_exists(file_path)

        file_size = os.path.getsize(file_path)
        logger.debug(f"File {file_path} has file size: {file_size}")
        return file_size

    @staticmethod
    def file_is_empty(file_path: str):
        file_size = File.get_file_size(file_path)
        return file_size == 0

    @staticmethod
    def copy_file(src_file_path: str, dst_file_path: str, overwrite: bool = True):
        """#TODO"""
        File.validate_file_exists(src_file_path)
        File.validate_path(dst_file_path)

        #
        if os.path.isfile(dst_file_path):
            if overwrite:
                os.remove(dst_file_path)
                shutil.copyfile(
                    src_file_path,
                    dst_file_path,
                )
                logger.debug(f"File {dst_file_path} overwritten by {src_file_path}.")
        else:
            logger.debug(f"File {src_file_path} copied to {dst_file_path}")
            shutil.copyfile(
                src_file_path,
                dst_file_path,
            )

    @staticmethod
    def delete_file(file_path: str):
        File.validate_path(file_path)
        if File.file_exists(file_path):
            os.remove(file_path)
            logger.debug(f"Deleted file {file_path}.")
        else:
            logger.debug(f"Tried to delete file {file_path} but it doesn't exist.")

    @staticmethod
    def files_differ(file_path_1: str, file_path_2: str, verbose: bool = False):
        """#TODO"""
        File.validate_file_exists(file_path_1)
        File.validate_file_exists(file_path_2)

        with open(file_path_1, "r") as f:
            a = set(f.readlines())
        with open(file_path_2, "r") as f:
            b = set(f.readlines())

        diff = [f"-\t{x}" if x in a else f"+\t{x}" for x in list(a ^ b)]

        if verbose:
            diff_str = "\n".join(diff)
            logger.debug(f"Diff between {file_path_1} and {file_path_2}:\n{diff_str}")

        return len(diff) != 0

    @staticmethod
    def file_exists(file_path: str):
        File.validate_path(file_path)
        return os.path.isfile(file_path)

    @staticmethod
    def read_file_lines(file_path: str):
        with open(file_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]
        return lines

    @staticmethod
    def file_is_gzipped(file_path: str):
        with open(file_path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"

    @staticmethod
    def validate_file_exists(file_path: str):
        if not File.file_exists(file_path):
            raise FileNotFoundError(f"File {file_path} does not exist.")

    @staticmethod
    def validate_file_is_not_empty(file_path: str):
        File.validate_file_exists(file_path)
        if File.file_is_empty(file_path):
            raise Exception(f"File {file_path} is empty.")

    def open_file(self) -> Union[TextIO, gzip.GzipFile, tarfile.TarFile]:
        _, file_extension = os.path.splitext(self.path)

        if file_extension == '.tar.gz' or file_extension == '.tgz':  # tarball
            return tarfile.open(self.path, "r:gz")
        elif file_extension == '.gz':  # gzip
            return gzip.open(self.path, 'rt')  # 'rt' mode for reading text files
        else:
            return open(self.path, 'r')


class SMIFile(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def read_dataframe(self):
        return self.read_dataframe_from_smi_file(self.path)

    @staticmethod
    def read_dataframe_from_smi_file(smi_file_path: str):
        File.validate_file_exists(smi_file_path)

        #
        data = []
        with open(smi_file_path, "r") as f:
            for line in f.readlines():
                line_elements = line.strip().split()
                if len(line_elements) != 2:
                    raise Exception(
                        f"Line in .smi file does not contain expected number of columns (2): {line_elements}"
                    )
                smiles_string, zinc_id = line_elements
                SMIFile.validate_smiles_string(smiles_string)
                # TODO: validate zinc_id
                data.append(
                    {
                        "zinc_id": zinc_id,
                        "smiles": smiles_string,
                    }
                )
        df = pd.DataFrame.from_records(data)

        return df

    @staticmethod
    def validate_smiles_string(smiles_string: str):
        m = Chem.MolFromSmiles(smiles_string, sanitize=False)
        if m is None:
            raise Exception(f"Invalid SMILES: {smiles_string}")
        else:
            try:
                Chem.SanitizeMol(m)
            except:
                raise Exception(f"Invalid chemistry in SMILES: {smiles_string}")


class SDIFile(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def write_tgz(self, tgz_file_name: str, archive_dir_name: Optional[str] = None, filter_regex: str = "(.*?)"):
        if archive_dir_name is None:
            archive_dir_name = File.get_file_name_of_file(tgz_file_name)
            archive_dir_name = re.sub(".tgz$", "", archive_dir_name)
            archive_dir_name = re.sub(".tar.gz$", "", archive_dir_name)
        db2_file_paths = self.read_lines()
        pattern = re.compile(filter_regex)
        temp_dir_name = str(uuid.uuid4())
        os.mkdir(temp_dir_name)
        with tarfile.open(tgz_file_name, "w:gz") as tar:
            i = 0
            for db2_file_path in db2_file_paths:
                if pattern.match(db2_file_path):
                    dst_file_name = f"{i+1}.db2"
                    dst_file_path = os.path.join(temp_dir_name, dst_file_name)
                    file_path_in_archive = os.path.join(archive_dir_name, dst_file_name)
                    if File.file_is_gzipped(db2_file_path):
                        with gzip.open(db2_file_path, "rb") as f_in:
                            with open(dst_file_path, "wb") as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    else:
                        with open(db2_file_path, "r") as f_in:
                            with open(dst_file_path, "w") as f_out:
                                shutil.copyfileobj(f_in, f_out)
                    tar.add(dst_file_path, arcname=file_path_in_archive)
                    i += 1
        shutil.rmtree(temp_dir_name)


class ProgramFile(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)


class LogFile(File):
    """#TODO"""

    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)


class IndockFile(File):
    """
    The INDOCK file is the main parameters file for the DOCK program.

    `indock_file_generation_dict` corresponds to the key-values of blastermaster_config.yaml
    """

    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def write(
        self,
        dock_files,
        indock_file_generation_dict,
        dock_files_dir_name="dockfiles",
        flex_groups=None,
        flex_0_file=None,
        use_flex=False,
        flexible_penalty_m=None,
    ):
        """takes a bunch of config, writes an appropriate INDOCK file"""
        if flex_groups is None:
            flex_groups = []

        def get_yes_or_no(boolean):
            if boolean:
                return "yes"
            else:
                return "no"

        #
        File.validate_file_exists(dock_files.electrostatics_phi_size_file.path)
        with open(dock_files.electrostatics_phi_size_file.path, "r") as f:
            try:
                phi_size = int(f.readline().strip())
            except Exception as e:
                raise Exception(
                    "Problem encountered while reading electrostatics phi size file. Check electrostatics phi size file."
                )

        # TODO: parametrize dock version
        header = f"""DOCK 3.8 parameter
#####################################################
# NOTE: split_database_index is reserved to specify a list of files
# defults for large scale docking.
ligand_atom_file               {indock_file_generation_dict['ligand_atom_file']}
#####################################################
#                             OUTPUT
output_file_prefix            {indock_file_generation_dict['output_file_prefix']}
#####################################################
#                             MATCHING
match_method                  {indock_file_generation_dict['match_method']}
distance_tolerance            {indock_file_generation_dict['distance_tolerance']}
match_goal                    {indock_file_generation_dict['match_goal']}
distance_step                 {indock_file_generation_dict['distance_step']}
distance_maximum              {indock_file_generation_dict['distance_maximum']}
timeout                       {indock_file_generation_dict['timeout']}
nodes_maximum                 {indock_file_generation_dict['nodes_maximum']}
nodes_minimum                 {indock_file_generation_dict['nodes_minimum']}
bump_maximum                  {indock_file_generation_dict['bump_maximum']}
bump_rigid                    {indock_file_generation_dict['bump_rigid']}
mol2_score_maximum            {indock_file_generation_dict['mol2_score_maximum']}
#####################################################
#                             COLORING
chemical_matching             {get_yes_or_no(indock_file_generation_dict['chemical_matching'])}
case_sensitive                {get_yes_or_no(indock_file_generation_dict['case_sensitive'])}
#####################################################
#                             SEARCH MODE
atom_minimum                  {indock_file_generation_dict['atom_minimum']}
atom_maximum                  {indock_file_generation_dict['atom_maximum']}
number_save                   {indock_file_generation_dict['number_save']}
number_write                  {indock_file_generation_dict['number_write']}
flush_int                     {indock_file_generation_dict['flush_int']}
#molecules_maximum            100000
check_clashes                 {get_yes_or_no(indock_file_generation_dict['check_clashes'])}
do_premax                     {get_yes_or_no(indock_file_generation_dict['do_premax'])}
do_clusters                   {get_yes_or_no(indock_file_generation_dict['do_clusters'])}
#####################################################
#                             SCORING
ligand_desolvation            {indock_file_generation_dict['ligand_desolvation']}
#vdw_maximum                   1.0e10
ligand_desolv_scale           {indock_file_generation_dict['ligand_desolv_scale']}
electrostatic_scale           {indock_file_generation_dict['electrostatic_scale']}
vdw_scale                     {indock_file_generation_dict['vdw_scale']}
internal_scale                {indock_file_generation_dict['internal_scale']}
per_atom_scores               {get_yes_or_no(indock_file_generation_dict['per_atom_scores'])}
##################################################### 
#                             DOCKovalent 
dockovalent                   {get_yes_or_no(indock_file_generation_dict['dockovalent'])}
bond_len                      {indock_file_generation_dict['bond_len']}
bond_ang1                     {indock_file_generation_dict['bond_ang1']}
bond_ang2                     {indock_file_generation_dict['bond_ang2']}
len_range                     {indock_file_generation_dict['len_range']}
len_step                      {indock_file_generation_dict['len_step']}
ang1_range                    {indock_file_generation_dict['ang1_range']}
ang2_range                    {indock_file_generation_dict['ang2_range']}
ang1_step                     {indock_file_generation_dict['ang1_step']}
ang2_step                     {indock_file_generation_dict['ang2_step']}
#####################################################
#                    MINIMIZATION
minimize                      {get_yes_or_no(indock_file_generation_dict['minimize'])}
sim_itmax                     {indock_file_generation_dict['sim_itmax']}
sim_trnstep                   {indock_file_generation_dict['sim_trnstep']}
sim_rotstep                   {indock_file_generation_dict['sim_rotstep']}
sim_need_to_restart           {indock_file_generation_dict['sim_need_to_restart']}
sim_cnvrge                    {indock_file_generation_dict['sim_cnvrge']}
min_cut                       {indock_file_generation_dict['min_cut']}
iseed                         {indock_file_generation_dict['iseed']}
##################################################### 
##                 Monte Carlo OPTIMIZATION
#monte_carlo                   no 
#mc_itmax                      500
#mc_accpt                      250
#mc_temp                       298.15
#mc_trnstep                    0.2
#mc_rotstep                    5.0
#mc_iseed                      777
#####################################################
# INPUT FILES / THINGS THAT CHANGE
"""

        #
        with open(self.path, "w") as f:
            f.write(header)
            f.write(
                f"receptor_sphere_file          {os.path.join('..', dock_files_dir_name, dock_files.matching_spheres_file.name)}\n"
            )
            f.write(
                f"vdw_parameter_file            {os.path.join('..', dock_files_dir_name, dock_files.vdw_parameters_file.name)}\n"
            )
            f.write(f"delphi_nsize                  {phi_size}\n")
            if not use_flex:  # normal docking, no flexible sidechains
                f.write(
                    f"flexible_receptor             {get_yes_or_no(indock_file_generation_dict['flexible_receptor'])}\n"
                )
                f.write(
                    f"total_receptors               {indock_file_generation_dict['total_receptors']}\n"
                )
                f.write("############## grids/data for one receptor\n")
                f.write(
                    f"rec_number                    {indock_file_generation_dict['rec_number']}\n"
                )
                f.write(
                    f"rec_group                     {indock_file_generation_dict['rec_group']}\n"
                )
                f.write(
                    f"rec_group_option              {indock_file_generation_dict['rec_group_option']}\n"
                )
                f.write(
                    f"solvmap_file                  {os.path.join('..', dock_files_dir_name, dock_files.ligand_desolvation_heavy_file.name)}\n"
                )
                f.write(
                    f"hydrogen_solvmap_file         {os.path.join('..', dock_files_dir_name, dock_files.ligand_desolvation_hydrogen_file.name)}\n"
                )
                f.write(
                    f"delphi_file                   {os.path.join('..', dock_files_dir_name, dock_files.electrostatics_trim_phi_file.name)}\n"
                )
                f.write(
                    f"chemgrid_file                 {os.path.join('..', dock_files_dir_name, dock_files.vdw_file.name)}\n"
                )
                f.write(
                    f"bumpmap_file                  {os.path.join('..', dock_files_dir_name, dock_files.vdw_bump_map_file.name)}\n"
                )
                f.write("#####################################################\n")
                f.write("#                             STRAIN\n")
                f.write(
                    f"check_strain                  {get_yes_or_no(indock_file_generation_dict['check_strain'])}\n"
                )
                f.write(
                    f"total_strain                  {indock_file_generation_dict['total_strain']}\n"
                )
                f.write(
                    f"max_strain                    {indock_file_generation_dict['max_strain']}\n"
                )
                f.write("############## end of INDOCK\n")
            else:  # flexible docking
                raise NotImplementedError
                # TODO
                """
                # flex_groups contains relevant data
                total_groups = [len(sub_group) for sub_group in flex_groups]
                f.write("flexible_receptor             yes\n")
                f.write("score_each_flex               yes\n")
                f.write(f"total_receptors               {str(total_groups)}\n")
                for i, sub_group in enumerate(flex_groups):
                    for j, one_rec in enumerate(sub_group):
                        energy = util.occupancy_to_energy(one_rec[1], flexible_penalty_m)
                        f.write("############## grids/data for one receptor\n")
                        f.write(f"## residues: {str(one_rec[2])}{one_rec[3]}\n")
                        f.write(f"## occupancy: {str(one_rec[1])}\n")
                        f.write(
                            f"## energy: {str(util.occupancy_to_energy(one_rec[1], 1.0))}\n"
                        )
                        f.write(f"## multiplier: {str(flexible_penalty_m)}\n")
                        f.write(f"## penalty: {str(energy)}\n")
                        f.write(f"rec_number                    {str(one_rec[0])}\n")
                        f.write(f"rec_group                     {str(i + 1)}\n")
                        f.write(f"rec_group_option              {str(j + 1)}\n")
                        f.write(f"rec_energy                    {str(energy)}\n")
                        f.write(
                            f"solvmap_file                  {os.path.join(output_dir.path, str(one_rec[0]), ligand_desolvation_pdb_outfile.path + ligand_desolvation_heavy_name)}\n"
                        )
                        f.write(
                            f"hydrogen_solvmap_file         {os.path.join(output_dir.path, str(one_rec[0]), ligand_desolvation_pdb_outfile.path + ligand_desolvation_hydrogen_name)}\n"
                        )
                        if (flex_0_file is None) or (
                            (str(one_rec[0]), "electrostatics")
                            not in flex_0_file  # TODO: this will cause an error since flex_0_file is potentially None
                        ):
                            f.write(
                                f"delphi_file                   {os.path.join(output_dir.path, str(one_rec[0]), dock_files.electrostatics_trim_phi_file.path)}\n"
                            )
                        else:  # means we are using implicitly 0 grid here, just write 0
                            f.write("delphi_file                   0\n")
                        f.write(
                            f"chemgrid_file                 {os.path.join(output_dir.path, str(one_rec[0]), vdw_prefix + '.vdw')}\n"
                        )
                        f.write(
                            f"bumpmap_file                  {os.path.join(output_dir.path, str(one_rec[0]), vdw_prefix + '.bmp')}\n"
                        )
                f.write("############## end of INDOCK\n")
                """


class TarballFile(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def iterate_over_tarball_member_files(self) -> Generator[tarfile.TarInfo, None, None]:
        with TarballFile(self.path).open_file() as tar:
            for member in tar.getmembers():
                if member.isfile():
                    yield member

    def extract(self, extraction_dir_path: str) -> None:
        if not self.exists:
            raise Exception(f"Tarball `{self.path}` does not exist.")

        if not Dir(extraction_dir_path).exists:
            raise Exception(f"Directory `{extraction_dir_path}` does not exist.")

        #
        with self.open_file() as tar:
            tar.extractall(path=extraction_dir_path)


class DB2File(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def get_molecule_name(self) -> Optional[str]:
        with self.open_file() as f:
            for line in f:
                columns = line.split()
                if columns[0] == "M":
                    return columns[1]
        return None


class OutdockFile(File):

    COLUMN_NAMES = [  # TODO: This depends on the version of DOCK 3 being used. Figure out how to make this more robust.
        "mol#",
        "id_num",
        "flexiblecode",
        "matched",
        "nscored",
        "time",
        "hac",
        "setnum",
        "matnum",
        "rank",
        "charge",
        "elect",
        "gist",
        "vdW",
        "psol",
        "asol",
        "tStrain",
        "mStrain",
        "rec_d",
        "r_hyd",
        "Total",
    ]

    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

    def get_dataframe(self):
        File.validate_file_exists(self.path)
        with open(self.path, "r", errors="ignore") as f:
            #
            lines = [x.strip() for x in f.readlines()]

            #
            if not lines[-1].startswith("elapsed time (sec):"):  # TODO: This depends on the version of DOCK 3 being used. Figure out how to make this more robust.
                raise Exception(f"Final line of OutdockFile {self.path} does not begin with 'elapsed time (sec):', indicating a failure of some kind.")

            # find first ligand line ("Input ligand: [...]")
            first_db2_line_index = None
            for i, line in enumerate(lines):
                if line.strip().endswith(".db2") or line.strip().endswith(".db2.gz"):
                    first_db2_line_index = i
                    break

            #
            header_line_index = None
            for i, line in enumerate(lines):
                if line.strip().startswith(self.COLUMN_NAMES[0]) and line.strip().endswith(self.COLUMN_NAMES[-1]):  # TODO: unfortunately this brittle solution will have to do for now, since this is the only way to be compatible with both DOCK 3.7 and 3.8
                    header_line_index = i
                    break
            if header_line_index is None:
                raise Exception(
                    f"Header line not found when reading OutdockFile: {self.path}"
                )

            #
            lines = [lines[first_db2_line_index]] + lines[header_line_index + 1 :]

            #
            db2_file_line_indices = [
                i
                for i, line in enumerate(lines)
                if line.endswith(".db2") or line.endswith(".db2.gz")
            ]  # TODO: this is quite brittle. find a better way.
            if len(db2_file_line_indices) % 2 != 0:
                raise Exception(f"Cannot parse OutdockFile: {self.path}")

            #
            open_file_line_indices = [
                db2_file_line_indices[i]
                for i in range(len(db2_file_line_indices))
                if i % 2 == 0
            ]
            close_file_line_indices = [
                db2_file_line_indices[i]
                for i in range(len(db2_file_line_indices))
                if i % 2 == 1
            ]

            #
            new_close_file_line_indices = []
            for close_file_line_index in close_file_line_indices:
                new_close_file_line_index = close_file_line_index
                if (
                    "close the file:" not in lines[close_file_line_index]
                ):  # apparently necessary if path is too long b/c Fortran is wacky
                    new_close_file_line_index -= 1
                new_close_file_line_indices.append(close_file_line_index)
            close_file_line_indices = new_close_file_line_indices

            #
            data = []
            db2_file_paths = []
            df_column_names = ["db2_file_path"] + self.COLUMN_NAMES
            for open_file_line_index, close_file_line_index in zip(
                open_file_line_indices, close_file_line_indices
            ):
                db2_file_path = (
                    lines[open_file_line_index]
                    .replace("open the file:", "")
                    .replace("Input ligand:", "")
                    .strip()
                )
                db2_file_paths.append(db2_file_path)
                for data_row_line_index in range(
                    open_file_line_index + 1, close_file_line_index
                ):
                    data_row_line = lines[data_row_line_index]
                    data_row = data_row_line.strip().split()
                    if data_row[0].isdigit():
                        data_row = [db2_file_path] + data_row
                    else:
                        data_row = [db2_file_path] + [
                            np.nan for _ in range(len(self.COLUMN_NAMES))
                        ]

                    # pad missing columns with NaN
                    if len(data_row) != len(df_column_names):
                        num_missing = len(df_column_names) - len(data_row)
                        data_row += [np.nan for _ in range(num_missing)]

                    #
                    data_row_dict = {
                        column_name: data_row[i]
                        for i, column_name in enumerate(df_column_names)
                    }
                    data.append(data_row_dict)

            return pd.DataFrame.from_records(data)


class Mol2Headers(Enum):
    ALT_TYPE = "@<TRIPOS>ALT_TYPE"
    ANCHOR_ATOM = "@<TRIPOS>ANCHOR_ATOM"
    ASSOCIATED_ANNOTATION = "@<TRIPOS>ASSOCIATED_ANNOTATION"
    ATOM = "@<TRIPOS>ATOM"
    BOND = "@<TRIPOS>BOND"
    CENTER_OF_MASS = "@<TRIPOS>CENTER_OF_MASS"
    CENTROID = "@<TRIPOS>CENTROID"
    COMMENT = "@<TRIPOS>COMMENT"
    CRYSIN = "@<TRIPOS>CRYSIN"
    DICT = "@<TRIPOS>DICT"
    DATA_FILE = "@<TRIPOS>DATA_FILE"
    EXTENSION_POINT = "@<TRIPOS>EXTENSION_POINT"
    FF_PBC = "@<TRIPOS>FF_PBC"
    FFCON_ANGLE = "@<TRIPOS>FFCON_ANGLE"
    FFCON_DIST = "@<TRIPOS>FFCON_DIST"
    FFCON_MULTI = "@<TRIPOS>FFCON_MULTI"
    FFCON_RANGE = "@<TRIPOS>FFCON_RANGE"
    FFCON_TORSION = "@<TRIPOS>FFCON_TORSION"
    LINE = "@<TRIPOS>LINE"
    LSPLANE = "@<TRIPOS>LSPLANE"
    MOLECULE = "@<TRIPOS>MOLECULE"
    NORMAL = "@<TRIPOS>NORMAL"
    QSAR_ALIGN_RULE = "@<TRIPOS>QSAR_ALIGN_RULE"
    RING_CLOSURE = "@<TRIPOS>RING_CLOSURE"
    ROTATABLE_BOND = "@<TRIPOS>ROTATABLE_BOND"
    SEARCH_DIST = "@<TRIPOS>SEARCH_DIST"
    SEARCH_OPTIONS = "@<TRIPOS>SEARCH_OPTIONS"
    SET = "@<TRIPOS>SET"
    SUBSTRUCTURE = "@<TRIPOS>SUBSTRUCTURE"
    U_FEAT = "@<TRIPOS>U_FEAT"
    UNITY_ATOM_ATTR = "@<TRIPOS>UNITY_ATOM_ATTR"
    UNITY_BOND_ATTR = "@<TRIPOS>UNITY_BOND_ATTR"


MOL2_HEADER_INDICATOR = "@<TRIPOS>"
MOL2_HEADER_STARTING_MOL2_BLOCK = Mol2Headers.MOLECULE.value


def find_nth_instance_of_line_starting_with_substring(lines: List[str], substring: str, n: int) -> Optional[int]:
    """
    Finds the index of the nth line that starts with the given substring.

    Parameters:
        lines (List[str]): List of lines to search.
        substring (str): Substring to find.
        n (int): The nth instance of the substring to find.

    Returns:
        Optional[int]: The index of the nth line that starts with the substring, or None if not found.
    """

    if n < 1:
        raise ValueError("n must be >= 1")

    num_found = 0
    for i, line in enumerate(lines):
        if line.startswith(substring):
            num_found += 1
            if num_found == n:
                return i
    return None


def remove_leading_invalid_mol2_lines(mol2_lines: List[str]) -> List[str]:
    """Removes leading lines that are invalid for a mol2 block from a list of mol2 lines."""

    new_lines = copy.copy(mol2_lines)
    while len(new_lines) > 0:
        first_line = new_lines[0].strip()
        if first_line.startswith(MOL2_HEADER_INDICATOR) or first_line.startswith(
                "#"):  # only valid start to a mol2 block
            break
        new_lines.pop(0)

    return new_lines


def get_leading_comment_block_end_index(lines: List[str]) -> Union[int, None]:
    """Get the end index of the leading comment block if it exists."""

    leading_comment_block_end_index = None
    found_leading_comment_block = False
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith("#"):
            found_leading_comment_block = True
            leading_comment_block_end_index = i
        elif line == "":
            if found_leading_comment_block:
                leading_comment_block_end_index = i
            continue
        elif line != "":
            break

    return leading_comment_block_end_index


def get_trailing_comment_block_start_index(lines: List[str]) -> Union[int, None]:
    """Get the start index of the trailing comment block if it exists."""

    trailing_comment_block_start_index = None
    for i in reversed(range(len(lines))):
        line = lines[i].strip()
        if line.startswith("#"):
            trailing_comment_block_start_index = i
        elif line == "":
            continue
        elif line != "":
            break

    return trailing_comment_block_start_index


def extract_leading_comment_block(lines: List[str]) -> List[str]:
    """Extracts leading comment lines from a list of lines."""

    leading_comment_end = get_leading_comment_block_end_index(lines)
    if leading_comment_end is not None:
        leading_comment_lines = lines[:leading_comment_end + 1]
    else:
        leading_comment_lines = []

    return leading_comment_lines


def remove_leading_comment_block(lines: List[str]) -> List[str]:
    """Removes the leading comment lines from a list of lines."""

    leading_comment_end = get_leading_comment_block_end_index(lines)
    if leading_comment_end is not None:
        del lines[:leading_comment_end + 1]

    return lines


def extract_trailing_comment_block(lines: List[str]) -> List[str]:
    """Extracts trailing comment lines from a list of lines."""

    trailing_comment_start = get_trailing_comment_block_start_index(lines)
    if trailing_comment_start is not None:
        trailing_comment_lines = lines[trailing_comment_start:]
    else:
        trailing_comment_lines = []

    return trailing_comment_lines


def remove_trailing_comment_block(lines: List[str]) -> List[str]:
    """Removes the trailing comment lines from a list of lines."""

    trailing_comment_start = get_trailing_comment_block_start_index(lines)
    if trailing_comment_start is not None:
        del lines[trailing_comment_start:]

    return lines


class Mol2DataRecord(object):
    def __init__(self, record_lines: List[str]):
        header, data_rows = self.parse_record_lines(record_lines)
        self.header = header
        self.data_rows = data_rows

    def __str__(self):
        return get_text_block(
            rows=self.data_rows,
            header=self.header,
            column_alignment='left',
            num_spaces_between_columns=4,
            num_spaces_before_line=4,
        )

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def parse_record_lines(record_lines: List[str]) -> Tuple[str, List[List[str]]]:
        header = record_lines[0].strip()
        data_rows = [line.split() for line in record_lines[1:]]
        return header, data_rows

    @staticmethod
    def line_is_valid_starting_line_to_mol2_data_record(line: str) -> bool:
        """Returns True if the line is a valid starting line to a mol2 data record."""
        return line.strip().startswith(MOL2_HEADER_INDICATOR) or line.strip().startswith("#")


class Mol2Block(object):
    def __init__(self, mol2_block_lines: List[str]) -> None:
        leading_comment_block_end_index = get_leading_comment_block_end_index(mol2_block_lines)
        if leading_comment_block_end_index is None:
            self.comment_lines = []
            remaining_lines = mol2_block_lines
        else:
            if len(mol2_block_lines) <= leading_comment_block_end_index + 1:
                raise ValueError(f"No records found in the supplied mol2 block lines: {mol2_block_lines}")
            self.comment_lines = mol2_block_lines[:leading_comment_block_end_index + 1]
            remaining_lines = mol2_block_lines[leading_comment_block_end_index + 1:]

        # Validate the records
        mol2_headers = [h.value for h in Mol2Headers]
        for header in mol2_headers:
            count = len([line for line in remaining_lines if line.strip() == header])
            if count > 1:
                raise ValueError(f"More than one data record of type {header} found in the supplied mol2 block lines.")

        # Parse the records
        records_dict = collections.OrderedDict()
        for record in self.split_mol2_block_into_data_records(remaining_lines):
            record_type_name = record.header.replace(MOL2_HEADER_INDICATOR, '').strip()  # e.g., "@<TRIPOS>ATOM" -> "ATOM"
            records_dict[record_type_name] = record
        self.data_records = records_dict

    def __str__(self):
        return "\n".join(self.comment_lines + [str(record) for record in self.data_records.values()])

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def split_mol2_block_into_data_records(mol2_block_lines: List[str]) -> List[Mol2DataRecord]:
        """Splits the supplied mol2 block into its component data records."""

        if not Mol2DataRecord.line_is_valid_starting_line_to_mol2_data_record(mol2_block_lines[0]):
            raise ValueError("The first line of the supplied mol2 block is not a valid record header or comment.")

        records = []
        n = 1
        should_break = False
        while True:
            start = find_nth_instance_of_line_starting_with_substring(mol2_block_lines, MOL2_HEADER_INDICATOR, n)
            end = find_nth_instance_of_line_starting_with_substring(mol2_block_lines, MOL2_HEADER_INDICATOR, n + 1)
            if end is None:
                end = len(mol2_block_lines)
                should_break = True

            record_lines = extract_trailing_comment_block(mol2_block_lines[:start]) + remove_trailing_comment_block(
                mol2_block_lines[start:end])
            records.append(Mol2DataRecord(record_lines))

            if should_break:
                break
            n += 1
        return records


class Mol2File(File):
    def __init__(self, path: str, validate_existence: bool = False):
        super().__init__(path=path, validate_existence=validate_existence)

        self.blocks = self.read_mol2_blocks(self.path)

    @staticmethod
    def read_mol2_blocks(mol2_file_path) -> List[Mol2Block]:
        """
        Reads and parses Mol2 blocks from a file specified by the `self.path` attribute.

        This function assumes that the file at `self.path` is a Mol2 file and processes it to
        extract information about molecules, atoms, and bonds. It captures this information into
        a list of Mol2Block objects, each representing a molecule and its associated details.

        Returns:
        --------
        list of Mol2Block
            A list of Mol2Block objects.

        Assumptions:
        ------------
        - The file is well-formatted according to the Mol2 standard.
        - The Mol2Block class is available and properly defined.

        Example Usage:
        --------------
        ```
        mol2_file = Mol2File(path="some_file.mol2")
        mol2_blocks = mol2_file.read_mol2_blocks()
        ```

        Notes:
        ------
        This function relies on `self.path` attribute to determine the file to read.
        """

        with open(mol2_file_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]

        mol2_blocks = Mol2File.split_mol2_file_lines_into_mol2_blocks(lines)

        return mol2_blocks

    @staticmethod
    def split_mol2_file_lines_into_mol2_blocks(mol2_lines: List[str]) -> List[Mol2Block]:
        """Split the supplied mol2 file lines into mol2 blocks."""

        # Remove leading invalid lines
        mol2_lines = remove_leading_invalid_mol2_lines(mol2_lines)

        #
        if not Mol2DataRecord.line_is_valid_starting_line_to_mol2_data_record(mol2_lines[0]):
            raise ValueError("The first line of the supplied mol2 block is not a valid record header or comment.")

        #
        mol2_blocks = []
        n = 1
        should_break = False
        while True:
            start = find_nth_instance_of_line_starting_with_substring(mol2_lines, MOL2_HEADER_STARTING_MOL2_BLOCK, n)
            end = find_nth_instance_of_line_starting_with_substring(mol2_lines, MOL2_HEADER_STARTING_MOL2_BLOCK, n + 1)
            if end is None:
                end = len(mol2_lines)
                should_break = True

            mol2_block_lines = extract_trailing_comment_block(mol2_lines[:start]) + remove_trailing_comment_block(mol2_lines[start:end])
            mol2_blocks.append(Mol2Block(mol2_block_lines))

            if should_break:
                break

            n += 1

        return mol2_blocks

    def write_mol2_file_with_molecules_cloned_and_transformed(
        self,
        rotation_matrix: np.ndarray,
        translation_vector: np.ndarray,
        write_path: str,
        num_applications: int = 1,
        bidirectional: bool = False,
    ) -> None:
        """
        Writes a new Mol2 file with molecules that are cloned and transformed based on
        given rotation matrix and translation vector parameters.

        Parameters:
        -----------
        rotation_matrix : numpy.ndarray
            A 3x3 rotation matrix to transform the molecule's coordinates.

        translation_vector : numpy.ndarray
            A 1x3 translation vector to transform the molecule's coordinates.

        write_path : str
            The path where the new Mol2 file will be written.

        num_applications : int, optional (default=1)
            The number of times the transformation should be applied to each molecule.

        bidirectional : bool, optional (default=False)
            If True, applies the transformation in both directions (forward and inverse).

        Examples:
        ---------
        # Given that `self` is an instance of a class that contains this method
        ```
        rotation_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        translation_vector = np.array([1, 1, 1])
        self.write_mol2_file_with_molecules_cloned_and_transformed(
            rotation_matrix,
            translation_vector,
            'output.mol2',
            num_applications=2,
            bidirectional=True
        )
        ```
        """

        # Validate rotation_matrix
        if not isinstance(rotation_matrix, np.ndarray) or rotation_matrix.shape != (3, 3):
            raise ValueError("`rotation_matrix` must be a 3x3 numpy.ndarray")

        # Validate translation_vector
        if not isinstance(translation_vector, np.ndarray) or translation_vector.shape != (3,):
            raise ValueError("`translation_vector` must be a 1x3 numpy.ndarray")

        # Validate write_path
        if not isinstance(write_path, str):
            raise ValueError("`write_path` must be a string")

        # Validate num_applications
        if not isinstance(num_applications, int) or num_applications < 1:
            raise ValueError("`num_applications` must be an integer greater than or equal to 1")

        # Validate bidirectional
        if not isinstance(bidirectional, bool):
            raise ValueError("`bidirectional` must be a boolean")

        #
        def transform(xyz: np.ndarray, rot_mat: np.ndarray, transl_vec: np.ndarray) -> np.ndarray:
            """Transforms the coordinates using the given rotation matrix and translation vector."""

            return np.dot(rot_mat, xyz) + transl_vec

        def get_inverse_transform(rot_mat: np.ndarray, transl_vec: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            """Calculates the inverse transformation."""

            a = np.array([[1.0, 0.0, 0.0, 0.0]])
            for i in range(3):
                a = np.concatenate(
                    (
                        a,
                        np.array(
                            [np.concatenate((np.array([transl_vec[i]]), rot_mat[i, :]))]
                        ),
                    ),
                    axis=0,
                )
            a_inv = np.linalg.inv(a)
            rot_mat_inv = a_inv[1:, 1:]
            transl_vec_inv = a_inv[1:, 0]

            return rot_mat_inv, transl_vec_inv

        #
        if bidirectional:
            rotation_matrix_inv, translation_vector_inv = get_inverse_transform(
                rotation_matrix, translation_vector
            )

        #
        with open(write_path, "w") as f:
            new_mol2_blocks = []
            for mol2_block in copy.deepcopy(self.blocks):
                #
                new_molecule_record_data_rows = []
                new_molecule_record_data_rows.append(mol2_block.data_records[Mol2Headers.MOLECULE.name].data_rows[0])
                molecule_row = mol2_block.data_records[Mol2Headers.MOLECULE.name].data_rows[1]
                if bidirectional:
                    multiplier = (2 * num_applications) + 1
                else:
                    multiplier = num_applications + 1
                new_molecule_row = [
                    int(molecule_row[0]) * multiplier,
                    int(molecule_row[1]) * multiplier,
                ] + molecule_row[2:]
                new_molecule_record_data_rows.append(new_molecule_row)
                if len(mol2_block.data_records[Mol2Headers.MOLECULE.name].data_rows) > 2:
                    new_molecule_record_data_rows += mol2_block.data_records[Mol2Headers.MOLECULE.name].data_rows[2:]

                #
                atom_element_to_id_nums_dict = collections.defaultdict(list)
                atom_names = [atom_row[1] for atom_row in mol2_block.data_records[Mol2Headers.ATOM.name].data_rows]
                for atom_name in atom_names:
                    element, id_num = [
                        token for token in re.split(r"(\d+)", atom_name) if token
                    ]
                    atom_element_to_id_nums_dict[element].append(int(id_num))

                #
                num_atoms = len(mol2_block.data_records[Mol2Headers.ATOM.name].data_rows)
                new_atom_record_data_rows = []
                for atom_row in mol2_block.data_records[Mol2Headers.ATOM.name].data_rows:
                    new_atom_record_data_rows.append(atom_row)

                def apply_to_atoms(rot_mat: np.ndarray, transl_vec: np.ndarray, n: int, num_app: int) -> None:
                    """
                    Applies transformations to atoms and appends the transformed atoms to `new_atom_record_data_rows`.

                    Parameters:
                    -----------
                    rot_mat : np.ndarray
                        A 3x3 rotation matrix to apply to each atom's coordinates.

                    transl_vec : np.ndarray
                        A 1x3 translation vector to apply to each atom's coordinates.

                    n : int
                        An integer multiplier for generating new atom IDs and names.

                    num_app : int
                        The number of times the transformation should be applied to each atom.

                    Assumptions:
                    ------------
                    - The function assumes that `new_atom_record_data_rows` and `mol2_block`
                      are defined in the outer scope.
                    - The function also assumes that `num_atoms` and `atom_element_to_id_nums_dict` are
                      defined and appropriately populated in the outer scope.
                    """

                    for atom_row in mol2_block.data_records[Mol2Headers.ATOM.name].data_rows:
                        atom_id = atom_row[0]
                        new_atom_id = f"{int(atom_id) + (n * num_atoms)}"
                        atom_name = atom_row[1]
                        element, id_num = [
                            token for token in re.split(r"(\d+)", atom_name) if token
                        ]
                        new_atom_name = f"{element}{int(id_num) + (n * max(atom_element_to_id_nums_dict[element]))}"
                        current_xyz = np.array(
                            [float(coord) for coord in atom_row[2:5]]
                        )
                        for j in range(num_app):
                            new_xyz = transform(current_xyz, rot_mat, transl_vec)
                            current_xyz = new_xyz
                        new_atom_row = (
                            [new_atom_id, new_atom_name] + list(new_xyz) + atom_row[5:]
                        )
                        new_atom_record_data_rows.append(new_atom_row)

                #
                for i, n in enumerate(list(range(1, num_applications + 1))):
                    apply_to_atoms(
                        rotation_matrix, translation_vector, n, num_app=i + 1
                    )

                #
                if bidirectional:
                    for i, n in enumerate(
                        list(range(num_applications + 1, (2 * num_applications) + 1))
                    ):
                        apply_to_atoms(
                            rotation_matrix_inv,
                            translation_vector_inv,
                            n,
                            num_app=i + 1,
                        )

                #
                num_bonds = len(mol2_block.data_records[Mol2Headers.BOND.name].data_rows)
                new_bond_record_data_rows = []
                for bond_row in mol2_block.data_records[Mol2Headers.BOND.name].data_rows:
                    new_bond_record_data_rows.append(bond_row)

                def apply_to_bonds(n):
                    """
                    Applies transformations to bonds and appends them to `new_bond_record_data_rows`.

                    This function takes an integer `n` and uses it to generate new bond IDs and associated atom IDs.
                    Transformed bond records are appended to the `new_bond_record_data_rows` list, assumed to be defined
                    in the outer scope.

                    Parameters:
                    -----------
                    n : int
                        An integer multiplier for generating new bond IDs and associated atom IDs.

                    Assumptions:
                    ------------
                    - Assumes `new_bond_record_data_rows` and `mol2_block` are defined in the outer scope.
                    - Assumes `num_atoms` and `num_bonds` are defined and appropriately populated in the outer scope.
                    """

                    for bond_row in mol2_block.data_records[Mol2Headers.BOND.name].data_rows:
                        new_bond_row = (
                            [int(bond_row[1]) + (n * num_bonds)]
                            + [int(num) + (n * num_atoms) for num in bond_row[1:3]]
                            + bond_row[3:]
                        )
                        new_bond_record_data_rows.append(new_bond_row)

                #
                for n in range(1, num_applications + 1):
                    apply_to_bonds(n)

                #
                if bidirectional:
                    for n in range(num_applications + 1, (2 * num_applications) + 1):
                        apply_to_bonds(n)

                #
                mol2_block.data_records[Mol2Headers.MOLECULE.name].data_rows = new_molecule_record_data_rows
                mol2_block.data_records[Mol2Headers.ATOM.name].data_rows = new_atom_record_data_rows
                mol2_block.data_records[Mol2Headers.BOND.name].data_rows = new_bond_record_data_rows

                #
                new_mol2_blocks.append(mol2_block)

            #
            f.write(f"{self.get_mol2_blocks_as_string(new_mol2_blocks)}\n")

    def __str__(self):
        """
        Returns string representation of the Mol2RecordSet object as would be written to a .mol2 file.
        """
        return Mol2File.get_mol2_blocks_as_string(self.blocks)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def get_mol2_blocks_as_string(mol2_blocks: List[Mol2Block]) -> str:
        """
        Returns a string representation of a list of Mol2Block objects as would be written to a .mol2 file.
        """
        return "\n".join([str(mol2_block) for mol2_block in mol2_blocks])


def get_text_block(
    rows: List[List[Any]],
    header: Optional[str] = None,
    column_alignment: str = "none",
    num_spaces_between_columns: int = 1,
    num_spaces_before_line: int = 0,
):
    """
    Generates a formatted text block from a list of rows.

    Parameters:
    -----------
    rows : list of list of any
        A list of rows, each of which is a list of items to be converted to strings.

    header : str, optional (default=None)
        An optional header text to be included at the top of the text block.

    column_alignment : str, optional (default="none")
        The text alignment for each column. Options are "left", "right", and "none".

    num_spaces_between_columns : int, optional (default=1)
        The number of spaces to insert between each column.

    num_spaces_before_line : int, optional (default=0)
        The number of spaces to insert before each line of text.

    Returns:
    --------
    str
        A formatted text block as a single string.

    Examples:
    ---------
    >>> get_text_block([[1, 20, 300], [4000, 50, 6]], header="Header", column_alignment="right")
    'Header
       1  20 300
    4000  50   6
    '
    """

    #
    if column_alignment not in ["left", "right", "none"]:
        raise ValueError(
            f"Invalid column_alignment value '{column_alignment}'. Valid options are 'left', 'right', and 'none'."
        )

    #
    formatted_lines = []
    if len(rows) > 0:
        rows = [[str(token) for token in row] for row in rows]

        max_row_size = max([len(row) for row in rows])
        columns = [
            [row[i] if i < len(row) else "" for row in rows] for i in range(max_row_size)
        ]
        column_max_token_length_list = [
            max([len(token) for token in column]) for column in columns
        ]
        spacing_between_columns = " " * num_spaces_between_columns
        formatted_lines = []
        for row in rows:
            formatted_tokens = []
            for i, token in enumerate(row):
                if column_alignment == "left":
                    formatted_token = token.ljust(column_max_token_length_list[i])
                elif column_alignment == "right":
                    formatted_token = token.rjust(column_max_token_length_list[i])
                elif column_alignment == "none":
                    formatted_token = token
                formatted_tokens.append(formatted_token)
            spacing_before_line = num_spaces_before_line * " "
            formatted_line = spacing_before_line + spacing_between_columns.join(
                formatted_tokens
            )
            formatted_lines.append(formatted_line)

    text_block = ""
    if header:
        text_block += f"{header}\n"
    text_block += "\n".join(formatted_lines)

    return text_block
