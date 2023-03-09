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

import numpy as np
import pandas as pd
from rdkit import Chem

from pydock3.util import validate_variable_type


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
INDOCK_FILE_NAME = "INDOCK"


class FileSystemEntity(object):
    """E.g., file, directory, symlink"""

    def __init__(self, path):
        self.path = path

    def __str__(self):
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

    def __init__(self, path, create=False, reset=False):
        super().__init__(path=path)

        #
        if create:
            self.create(reset=reset)

    @property
    def name(self):
        return Dir.extract_dir_name_from_dir_path(self.path)

    @staticmethod
    def extract_dir_name_from_dir_path(dir_path):
        return os.path.basename(dir_path)

    @property
    def exists(self):
        return Dir.dir_exists(self.path)

    @property
    def validate_existence(self):
        if not self.exists:
            raise Exception(f"Dir {self.path} does not exist.")

    @staticmethod
    def dir_exists(dir_path):
        return os.path.isdir(dir_path)

    def create(self, reset=False):
        """#TODO"""
        if reset:
            self.delete()
            pathlib.Path(self.path).mkdir(parents=True)
            logger.info(f"Reset directory {self}.")
        else:
            if os.path.exists(self.path):
                logger.debug(
                    f"Tried to create directory {self} with reset=False but directory already exists."
                )
            else:
                pathlib.Path(self.path).mkdir(parents=True)
                logger.info(f"Created directory {self}")

    def delete(self):
        if os.path.exists(self.path):
            shutil.rmtree(self.path, ignore_errors=True)
            while os.path.isdir(self.path):  # TODO: hmmm. is this valid?
                pass
            logger.info(f"Deleted directory {self}.")

    def reset(self):
        self.create(reset=True)

    def copy_in_file(self, src_file_path, dst_file_name=None, overwrite=True):
        """#TODO"""
        File.validate_file_exists(src_file_path)

        if dst_file_name is None:
            dst_file_name = File.get_file_name_of_file(src_file_path)
        dst_file = File(path=os.path.join(self.path, dst_file_name))
        dst_file.copy_from(src_file_path=src_file_path, overwrite=overwrite)

        return dst_file

    @staticmethod
    def validate_obj_is_dir(obj):
        validate_variable_type(obj, allowed_instance_types=(Dir,))


class File(FileSystemEntity):
    """#TODO"""

    def __init__(self, path):
        super().__init__(path=path)

    @property
    def name(self):
        return File.get_file_name_of_file(self.path)

    @property
    def datetime_last_modified(self):
        return self.get_datetime_file_was_last_modified(self.path)

    @staticmethod
    def get_file_name_of_file(file_path):
        File.validate_path(file_path)
        return os.path.basename(os.path.abspath(file_path))

    @staticmethod
    def get_dir_path_of_file(file_path):
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

    def copy_from(self, src_file_path, overwrite=True):
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
    def get_datetime_file_was_last_modified(file_path):
        File.validate_file_exists(file_path)

        datetime_last_modified = datetime.fromtimestamp(os.stat(file_path).st_mtime)
        logger.debug(f"File {file_path} was last modified at: {datetime_last_modified}")
        return datetime_last_modified

    @staticmethod
    def get_file_size(file_path):
        File.validate_file_exists(file_path)

        file_size = os.path.getsize(file_path)
        logger.debug(f"File {file_path} has file size: {file_size}")
        return file_size

    @staticmethod
    def file_is_empty(file_path):
        file_size = File.get_file_size(file_path)
        return file_size == 0

    @staticmethod
    def copy_file(src_file_path, dst_file_path, overwrite=True):
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
    def delete_file(file_path):
        File.validate_path(file_path)
        if File.file_exists(file_path):
            os.remove(file_path)
            logger.debug(f"Deleted file {file_path}.")
        else:
            logger.debug(f"Tried to delete file {file_path} but it doesn't exist.")

    @staticmethod
    def files_differ(file_path_1, file_path_2, verbose=False):
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
    def file_exists(file_path):
        File.validate_path(file_path)
        return os.path.isfile(file_path)

    @staticmethod
    def read_file_lines(file_path):
        with open(file_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]
        return lines

    @staticmethod
    def file_is_gzipped(file_path):
        with open(file_path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"

    @staticmethod
    def validate_file_exists(file_path):
        if not File.file_exists(file_path):
            raise FileNotFoundError(f"File {file_path} does not exist.")

    @staticmethod
    def validate_file_is_not_empty(file_path):
        File.validate_file_exists(file_path)
        if File.file_is_empty(file_path):
            raise Exception(f"File {file_path} is empty.")


class SMIFile(File):
    def __init__(self, path):
        super().__init__(path=path)

    def read_dataframe(self):
        self.read_dataframe_from_smi_file(self.path)

    @staticmethod
    def read_dataframe_from_smi_file(smi_file_path):
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
    def validate_smiles_string(smiles_string):
        m = Chem.MolFromSmiles(smiles_string, sanitize=False)
        if m is None:
            raise Exception(f"Invalid SMILES: {smiles_string}")
        else:
            try:
                Chem.SanitizeMol(m)
            except:
                raise Exception(f"Invalid chemistry in SMILES: {smiles_string}")


class SDIFile(File):
    def __init__(self, path):
        super().__init__(path=path)

    def write_tgz(self, tgz_file_name, archive_dir_name=None, filter_regex="(.*?)"):
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
    def __init__(self, path):
        super().__init__(path=path)


class LogFile(File):
    """#TODO"""

    def __init__(self, path):
        super().__init__(path=path)


class IndockFile(File):
    """
    The INDOCK file is the main parameters file for the DOCK program.

    `config_param_dict` corresponds to the key-values of blastermaster_config.yaml
    """

    def __init__(self, path):
        super().__init__(path=path)

    def write(
        self,
        dock_files,
        config_param_dict,
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
ligand_atom_file               {config_param_dict['indock_file_generation.ligand_atom_file']}
#####################################################
#                             OUTPUT
output_file_prefix            {config_param_dict['indock_file_generation.output_file_prefix']}
#####################################################
#                             MATCHING
match_method                  {config_param_dict['indock_file_generation.match_method']}
distance_tolerance            {config_param_dict['indock_file_generation.distance_tolerance']}
match_goal                    {config_param_dict['indock_file_generation.match_goal']}
distance_step                 {config_param_dict['indock_file_generation.distance_step']}
distance_maximum              {config_param_dict['indock_file_generation.distance_maximum']}
timeout                       {config_param_dict['indock_file_generation.timeout']}
nodes_maximum                 {config_param_dict['indock_file_generation.nodes_maximum']}
nodes_minimum                 {config_param_dict['indock_file_generation.nodes_minimum']}
bump_maximum                  {config_param_dict['indock_file_generation.bump_maximum']}
bump_rigid                    {config_param_dict['indock_file_generation.bump_rigid']}
mol2_score_maximum            {config_param_dict['indock_file_generation.mol2_score_maximum']}
#####################################################
#                             COLORING
chemical_matching             {get_yes_or_no(config_param_dict['indock_file_generation.chemical_matching'])}
case_sensitive                {get_yes_or_no(config_param_dict['indock_file_generation.case_sensitive'])}
#####################################################
#                             SEARCH MODE
atom_minimum                  {config_param_dict['indock_file_generation.atom_minimum']}
atom_maximum                  {config_param_dict['indock_file_generation.atom_maximum']}
number_save                   {config_param_dict['indock_file_generation.number_save']}
number_write                  {config_param_dict['indock_file_generation.number_write']}
flush_int                     {config_param_dict['indock_file_generation.flush_int']}
#molecules_maximum            100000
check_clashes                 {get_yes_or_no(config_param_dict['indock_file_generation.check_clashes'])}
do_premax                     {get_yes_or_no(config_param_dict['indock_file_generation.do_premax'])}
do_clusters                   {get_yes_or_no(config_param_dict['indock_file_generation.do_clusters'])}
#####################################################
#                             SCORING
ligand_desolvation            {config_param_dict['indock_file_generation.ligand_desolvation']}
#vdw_maximum                   1.0e10
ligand_desolv_scale           {config_param_dict['indock_file_generation.ligand_desolv_scale']}
electrostatic_scale           {config_param_dict['indock_file_generation.electrostatic_scale']}
vdw_scale                     {config_param_dict['indock_file_generation.vdw_scale']}
internal_scale                {config_param_dict['indock_file_generation.internal_scale']}
per_atom_scores               {get_yes_or_no(config_param_dict['indock_file_generation.per_atom_scores'])}
##################################################### 
#                             DOCKovalent 
dockovalent                   {get_yes_or_no(config_param_dict['indock_file_generation.dockovalent'])}
bond_len                      {config_param_dict['indock_file_generation.bond_len']}
bond_ang1                     {config_param_dict['indock_file_generation.bond_ang1']}
bond_ang2                     {config_param_dict['indock_file_generation.bond_ang2']}
len_range                     {config_param_dict['indock_file_generation.len_range']}
len_step                      {config_param_dict['indock_file_generation.len_step']}
ang1_range                    {config_param_dict['indock_file_generation.ang1_range']}
ang2_range                    {config_param_dict['indock_file_generation.ang2_range']}
ang1_step                     {config_param_dict['indock_file_generation.ang1_step']}
ang2_step                     {config_param_dict['indock_file_generation.ang2_step']}
#####################################################
#                    MINIMIZATION
minimize                      {get_yes_or_no(config_param_dict['indock_file_generation.minimize'])}
sim_itmax                     {config_param_dict['indock_file_generation.sim_itmax']}
sim_trnstep                   {config_param_dict['indock_file_generation.sim_trnstep']}
sim_rotstep                   {config_param_dict['indock_file_generation.sim_rotstep']}
sim_need_to_restart           {config_param_dict['indock_file_generation.sim_need_to_restart']}
sim_cnvrge                    {config_param_dict['indock_file_generation.sim_cnvrge']}
min_cut                       {config_param_dict['indock_file_generation.min_cut']}
iseed                         {config_param_dict['indock_file_generation.iseed']}
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
                    f"flexible_receptor             {get_yes_or_no(config_param_dict['indock_file_generation.flexible_receptor'])}\n"
                )
                f.write(
                    f"total_receptors               {config_param_dict['indock_file_generation.total_receptors']}\n"
                )
                f.write("############## grids/data for one receptor\n")
                f.write(
                    f"rec_number                    {config_param_dict['indock_file_generation.rec_number']}\n"
                )
                f.write(
                    f"rec_group                     {config_param_dict['indock_file_generation.rec_group']}\n"
                )
                f.write(
                    f"rec_group_option              {config_param_dict['indock_file_generation.rec_group_option']}\n"
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
                    f"check_strain                  {get_yes_or_no(config_param_dict['indock_file_generation.check_strain'])}\n"
                )
                f.write(
                    f"total_strain                  {config_param_dict['indock_file_generation.total_strain']}\n"
                )
                f.write(
                    f"max_strain                    {config_param_dict['indock_file_generation.max_strain']}\n"
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


class OutdockFile(File):

    COLUMN_NAMES = [
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

    def __init__(self, path):
        super().__init__(path=path)

    def get_dataframe(self):
        File.validate_file_exists(self.path)
        with open(self.path, "r", errors="ignore") as f:
            #
            lines = [x.strip() for x in f.readlines()]

            #
            if not lines[-1].startswith("elapsed time (sec):"):
                raise Exception("Final line of OutdockFile does not begin with 'elapsed time (sec):', indicating a failure of some kind.")

            # find first ligand line ("Input ligand: [...]")
            first_db2_line_index = None
            for i, line in enumerate(lines):
                if line.strip().endswith(".db2") or line.strip().endswith(".db2.gz"):
                    first_db2_line_index = i
                    break

            #
            header_line_index = None
            for i, line in enumerate(lines):
                if all([column_name in line for column_name in self.COLUMN_NAMES]):
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


class Mol2Record(object):
    MOLECULE_RECORD_HEADER = "@<TRIPOS>MOLECULE"
    ATOM_RECORD_HEADER = "@<TRIPOS>ATOM"
    BOND_RECORD_HEADER = "@<TRIPOS>BOND"

    def __init__(
        self, comment_lines, molecule_record_lines, atom_record_lines, bond_record_lines
    ):
        self.comment_lines = comment_lines
        self.molecule_record_lines = molecule_record_lines
        self.atom_record_lines = atom_record_lines
        self.bond_record_lines = bond_record_lines


class Mol2File(File):
    def __init__(self, path):
        super().__init__(path=path)

    def read_mol2_records(self):

        with open(self.path, "r") as f:
            remaining_lines = [line.strip() for line in f.readlines()]

        mol2_records = []

        start_index = remaining_lines.index(Mol2Record.MOLECULE_RECORD_HEADER)
        while True:
            pre_start_comment_lines = [
                line for line in remaining_lines[:start_index] if line.startswith("#")
            ]

            try:
                end_index = (
                    remaining_lines[start_index + 1 :].index(
                        Mol2Record.MOLECULE_RECORD_HEADER
                    )
                    + 1
                )
                mol2_record_lines = remaining_lines[start_index:end_index]
                mol2_record_string = "\n".join(mol2_record_lines)
                molecule_record_string, remaining_string = mol2_record_string.replace(
                    Mol2Record.MOLECULE_RECORD_HEADER, ""
                ).split(Mol2Record.ATOM_RECORD_HEADER)
                atom_record_string, bond_record_string = remaining_string.split(
                    Mol2Record.BOND_RECORD_HEADER
                )
                molecule_record_lines = [
                    line.split() for line in molecule_record_string.split("\n") if line
                ]
                atom_record_lines = [
                    line.split() for line in atom_record_string.split("\n") if line
                ]
                bond_record_lines = [
                    line.split() for line in bond_record_string.split("\n") if line
                ]
                mol2_record = Mol2Record(
                    comment_lines=pre_start_comment_lines,
                    molecule_record_lines=molecule_record_lines,
                    atom_record_lines=atom_record_lines,
                    bond_record_lines=bond_record_lines,
                )
                mol2_records.append(mol2_record)
                remaining_lines = remaining_lines[end_index:]
            except ValueError:
                mol2_record_lines = remaining_lines[start_index:]
                mol2_record_string = "\n".join(mol2_record_lines)
                molecule_record_string, remaining_string = mol2_record_string.replace(
                    Mol2Record.MOLECULE_RECORD_HEADER, ""
                ).split(Mol2Record.ATOM_RECORD_HEADER)
                atom_record_string, bond_record_string = remaining_string.split(
                    Mol2Record.BOND_RECORD_HEADER
                )
                molecule_record_lines = [
                    line.split() for line in molecule_record_string.split("\n") if line
                ]
                atom_record_lines = [
                    line.split() for line in atom_record_string.split("\n") if line
                ]
                bond_record_lines = [
                    line.split() for line in bond_record_string.split("\n") if line
                ]
                mol2_record = Mol2Record(
                    comment_lines=pre_start_comment_lines,
                    molecule_record_lines=molecule_record_lines,
                    atom_record_lines=atom_record_lines,
                    bond_record_lines=bond_record_lines,
                )
                mol2_records.append(mol2_record)
                break

        return mol2_records

    def write_mol2_file_with_molecules_cloned_and_transformed(
        self,
        rotation_matrix,
        translation_vector,
        write_path,
        num_applications=1,
        bidirectional=False,
    ):

        #
        def transform(xyz, rot_mat, transl_vec):
            return np.dot(rot_mat, xyz) + transl_vec

        def get_inverse_transform(rot_mat, transl_vec):
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
        def get_record_text_block(
            rows,
            header,
            alignment="right",
            num_spaces_before_line=5,
            num_spaces_between_columns=2,
        ):
            return get_text_block(
                rows,
                header,
                alignment=alignment,
                num_spaces_before_line=num_spaces_before_line,
                num_spaces_between_columns=num_spaces_between_columns,
            )

        #
        mol2_records = self.read_mol2_records()

        #
        if bidirectional:
            rotation_matrix_inv, translation_vector_inv = get_inverse_transform(
                rotation_matrix, translation_vector
            )

        #
        with open(write_path, "w") as f:

            for mol2_record in mol2_records:
                #
                new_molecule_record_lines = []
                new_molecule_record_lines.append(mol2_record.molecule_record_lines[0])
                molecule_row = mol2_record.molecule_record_lines[1]
                if bidirectional:
                    multiplier = (2 * num_applications) + 1
                else:
                    multiplier = num_applications + 1
                new_molecule_row = [
                    int(molecule_row[0]) * multiplier,
                    int(molecule_row[1]) * multiplier,
                ] + molecule_row[2:]
                new_molecule_record_lines.append(new_molecule_row)

                #
                atom_element_to_id_nums_dict = collections.defaultdict(list)
                atom_names = [atom_row[1] for atom_row in mol2_record.atom_record_lines]
                for atom_name in atom_names:
                    element, id_num = [
                        token for token in re.split(r"(\d+)", atom_name) if token
                    ]
                    atom_element_to_id_nums_dict[element].append(int(id_num))

                #
                num_atoms = len(mol2_record.atom_record_lines)
                new_atom_record_lines = []
                for atom_row in mol2_record.atom_record_lines:
                    new_atom_record_lines.append(atom_row)

                def apply_to_atoms(rot_mat, transl_vec, n, num_app):
                    for atom_row in mol2_record.atom_record_lines:
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
                        new_atom_record_lines.append(new_atom_row)

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
                num_bonds = len(mol2_record.bond_record_lines)
                new_bond_record_lines = []
                for bond_row in mol2_record.bond_record_lines:
                    new_bond_record_lines.append(bond_row)

                def apply_to_bonds(n):
                    for bond_row in mol2_record.bond_record_lines:
                        new_bond_row = (
                            [int(bond_row[1]) + (n * num_bonds)]
                            + [int(num) + (n * num_atoms) for num in bond_row[1:3]]
                            + bond_row[3:]
                        )
                        new_bond_record_lines.append(new_bond_row)

                #
                for n in range(1, num_applications + 1):
                    apply_to_bonds(n)

                #
                if bidirectional:
                    for n in range(num_applications + 1, (2 * num_applications) + 1):
                        apply_to_bonds(n)

                #
                f.write("\n".join(mol2_record.comment_lines) + "\n")
                f.write(
                    get_record_text_block(
                        new_molecule_record_lines, Mol2Record.MOLECULE_RECORD_HEADER
                    )
                )
                f.write(
                    get_record_text_block(
                        new_atom_record_lines, Mol2Record.ATOM_RECORD_HEADER
                    )
                )
                f.write(
                    get_record_text_block(
                        new_bond_record_lines, Mol2Record.BOND_RECORD_HEADER
                    )
                )


def get_text_block(
    rows,
    header=None,
    alignment="left",
    num_spaces_between_columns=1,
    num_spaces_before_line=0,
):
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
            if alignment == "left":
                formatted_token = token.ljust(column_max_token_length_list[i])
            elif alignment == "right":
                formatted_token = token.rjust(column_max_token_length_list[i])
            else:
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
    text_block += "\n".join(formatted_lines) + "\n"

    return text_block
