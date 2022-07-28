import logging
import os
import shutil
import pathlib
from datetime import datetime

import numpy as np
import pandas as pd
import rdkit

from ucsfdock.util import validate_variable_type


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
            if os.path.exists(self.path):
                shutil.rmtree(self.path)
                while os.path.isdir(self.path):
                    pass
                pathlib.Path(self.path).mkdir(parents=True)
                logger.info(f"Reset directory {self}.")
            else:
                pathlib.Path(self.path).mkdir(parents=True)
                logger.info(f"Created directory {self}.")
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
            shutil.rmtree(self.path)
            while os.path.isdir(self.path):
                pass
            logger.info(f"Deleted directory {self}.")

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
        self.copy_file(src_file_path=src_file_path, dst_file_path=self.path, overwrite=overwrite)

    def delete(self):
        self.delete_file(self.path)

    def is_newer_than_file(self, file_path):
        File.validate_path(file_path)
        return self.get_datetime_file_was_last_modified(self.path) > self.get_datetime_file_was_last_modified(file_path)

    def differs_from_file(self, file_path, verbose=False):
        return self.files_differ(self.path, file_path, verbose=verbose)

    def validate_is_not_empty(self):
        if self.is_empty:
            raise Exception(f"File is empty: {self}")

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

        with open(file_path_1, 'r') as f:
            a = set(f.readlines())
        with open(file_path_2, 'r') as f:
            b = set(f.readlines())

        diff = [f'-\t{x}' if x in a else f'+\t{x}' for x in list(a ^ b)]

        if verbose:
            diff_str = '\n'.join(diff)
            logger.debug(f"Diff between {file_path_1} and {file_path_2}:\n{diff_str}")

        return len(diff) != 0

    @staticmethod
    def file_exists(file_path):
        File.validate_path(file_path)
        return os.path.isfile(file_path)

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
        with open(smi_file_path, 'r') as f:
            for line in f.readlines():
                line_elements = line.strip().split()
                if len(line_elements) != 2:
                    raise Exception(f"Line in .smi file does not contain expected number of columns (2): {line_elements}")
                smiles_string, zinc_id = line_elements
                SMIFile.validate_smiles_string(smiles_string)
                # TODO: validate zinc_id
                data.append({
                    'zinc_id': zinc_id,
                    'smiles': smiles_string,
                })
        df = pd.DataFrame.from_records(data)

        return df

    @staticmethod
    def validate_smiles_string(smiles_string):
        m = rdkit.Chem.MolFromSmiles(smiles_string, sanitize=False)
        if m is None:
            raise Exception(f"Invalid SMILES: {smiles_string}")
        else:
            try:
                rdkit.Chem.SanitizeMol(m)
            except:
                raise Exception(f"Invalid chemistry in SMILES: {smiles_string}")


class SDIFile(File):
    def __init__(self, path):
        super().__init__(path=path)


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
        flex_groups=None,
        flex_0_file=None,
        use_flex=False,
        flexible_penalty_m=None,
    ):
        """takes a bunch of config, writes an appropriate INDOCK file"""
        if flex_groups is None:
            flex_groups = []

        #
        File.validate_file_exists(dock_files.electrostatics_phi_size_file.path)
        with open(dock_files.electrostatics_phi_size_file.path, 'r') as f:
            try:
                phi_size = int(f.readline().strip())
            except Exception as e:
                raise Exception("Problem encountered while reading electrostatics phi size file. Check electrostatics phi size file.")

        # TODO: parametrize dock version
        header = f"""DOCK 3.8 parameter
#####################################################
# NOTE: split_database_index is reserved to specify a list of files
# defults for large scale docking.
ligand_atom_file               {config_param_dict['indock.ligand_atom_file']}
#####################################################
#                             OUTPUT
output_file_prefix            {config_param_dict['indock.output_file_prefix']}
#####################################################
#                             MATCHING
match_method                  {config_param_dict['indock.match_method']}
distance_tolerance            {config_param_dict['indock.distance_tolerance']}
match_goal                    {config_param_dict['indock.match_goal']}
distance_step                 {config_param_dict['indock.distance_step']}
distance_maximum              {config_param_dict['indock.distance_maximum']}
timeout                       {config_param_dict['indock.timeout']}
nodes_maximum                 {config_param_dict['indock.nodes_maximum']}
nodes_minimum                 {config_param_dict['indock.nodes_minimum']}
bump_maximum                  {config_param_dict['indock.bump_maximum']}
bump_rigid                    {config_param_dict['indock.bump_rigid']}
mol2_score_maximum            {config_param_dict['indock.mol2_score_maximum']}
#####################################################
#                             COLORING
chemical_matching             {config_param_dict['indock.chemical_matching']}
case_sensitive                {config_param_dict['indock.case_sensitive']}
#####################################################
#                             SEARCH MODE
atom_minimum                  {config_param_dict['indock.atom_minimum']}
atom_maximum                  {config_param_dict['indock.atom_maximum']}
number_save                   {config_param_dict['indock.number_save']}
number_write                  {config_param_dict['indock.number_write']}
flush_int                     {config_param_dict['indock.flush_int']}
#molecules_maximum            100000
check_clashes                 {config_param_dict['indock.check_clashes']}
do_premax                     {config_param_dict['indock.do_premax']}
do_clusters                   {config_param_dict['indock.do_clusters']}
#####################################################
#                             SCORING
ligand_desolvation            {config_param_dict['indock.ligand_desolvation']}
#vdw_maximum                   1.0e10
ligand_desolv_scale           {config_param_dict['indock.ligand_desolv_scale']}
electrostatic_scale           {config_param_dict['indock.electrostatic_scale']}
vdw_scale                     {config_param_dict['indock.vdw_scale']}
internal_scale                {config_param_dict['indock.internal_scale']}
per_atom_scores               {config_param_dict['indock.per_atom_scores']}
##################################################### 
#                             DOCKovalent 
dockovalent                   {config_param_dict['indock.dockovalent']}
bond_len                      {config_param_dict['indock.bond_len']}
bond_ang1                     {config_param_dict['indock.bond_ang1']}
bond_ang2                     {config_param_dict['indock.bond_ang2']}
len_range                     {config_param_dict['indock.len_range']}
len_step                      {config_param_dict['indock.len_step']}
ang1_range                    {config_param_dict['indock.ang1_range']}
ang2_range                    {config_param_dict['indock.ang2_range']}
ang1_step                     {config_param_dict['indock.ang1_step']}
ang2_step                     {config_param_dict['indock.ang2_step']}
#####################################################
#                    MINIMIZATION
minimize                      {config_param_dict['indock.minimize']}
sim_itmax                     {config_param_dict['indock.sim_itmax']}
sim_trnstep                   {config_param_dict['indock.sim_trnstep']}
sim_rotstep                   {config_param_dict['indock.sim_rotstep']}
sim_need_to_restart           {config_param_dict['indock.sim_need_to_restart']}
sim_cnvrge                    {config_param_dict['indock.sim_cnvrge']}
min_cut                       {config_param_dict['indock.min_cut']}
iseed                         {config_param_dict['indock.iseed']}
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
                f"receptor_sphere_file          {dock_files.matching_spheres_file.path}\n"
            )
            f.write(
                f"vdw_parameter_file            {dock_files.vdw_parameters_file.path}\n"
            )
            f.write(f"delphi_nsize                  {phi_size}\n")
            if not use_flex:  # normal docking, no flexible sidechains
                f.write(f"flexible_receptor             {config_param_dict['indock.flexible_receptor']}\n")
                f.write(f"total_receptors               {config_param_dict['indock.total_receptors']}\n")
                f.write("############## grids/data for one receptor\n")
                f.write(f"rec_number                    {config_param_dict['indock.rec_number']}\n")
                f.write(f"rec_group                     {config_param_dict['indock.rec_group']}\n")
                f.write(f"rec_group_option              {config_param_dict['indock.rec_group_option']}\n")
                f.write(
                    f"solvmap_file                  {dock_files.ligand_desolvation_heavy_file.path}\n"
                )
                f.write(
                    f"hydrogen_solvmap_file         {dock_files.ligand_desolvation_hydrogen_file.path}\n"
                )
                f.write(
                    f"delphi_file                   {dock_files.electrostatics_trim_phi_file.path}\n"
                )
                f.write(f"chemgrid_file                 {dock_files.vdw_file.path}\n")
                f.write(
                    f"bumpmap_file                  {dock_files.vdw_bump_map_file.path}\n"
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
        with open(self.path, 'r', errors='ignore') as f:
            #
            lines = [x.strip() for x in f.readlines()]

            # find first ligand line ("Input ligand: [...]")
            first_lig_line_index = None
            for i, line in enumerate(lines):
                if line.strip().startswith("Input ligand"):
                    first_lig_line_index = i
                    break

            #
            header_line_index = None
            for i, line in enumerate(lines):
                if all([column_name in line for column_name in self.COLUMN_NAMES]):
                    header_line_index = i
                    break
            if header_line_index is None:
                raise Exception(f"Header line not found when reading OutdockFile: {self.path}")

            #
            lines = [lines[first_lig_line_index]] + lines[header_line_index+1:]

            #
            open_file_line_indices = [i for i, line in enumerate(lines) if line.startswith("open the file:") or line.startswith("Input ligand:")]

            #
            close_file_line_indices = []
            new_open_file_line_indices = []
            for i in open_file_line_indices:
                open_file_line = lines[i]
                db2_file_path = open_file_line.replace("open the file:", "").replace("Input ligand:", "").strip()
                close_file_line_index = None
                for j, line in enumerate(lines[i:]):
                    if line.startswith("close the file:"):
                        close_file_line_index = i + j
                        if line.replace("close the file:", "").strip() != db2_file_path:
                            raise Exception(f"Open file line {i+1} and close file line {close_file_line_index+1} do not match in OutdockFile: {self.path}")
                        break
                if close_file_line_index is None:
                    raise Exception(f"Corresponding close file line not found for open file line {i+1} in OutdockFile {self.path} : {open_file_line}")
                new_open_file_line_indices.append(i)
                close_file_line_indices.append(close_file_line_index)
            open_file_line_indices = new_open_file_line_indices

            #
            if len(open_file_line_indices) != len(close_file_line_indices):
                raise Exception(f"# of open file lines and # of close file lines do not match in OutdockFile: {self.path}")

            #
            db2_file_paths = []
            data = []
            df_column_names = ["db2_file_path"] + self.COLUMN_NAMES
            for open_file_line_index, close_file_line_index in zip(open_file_line_indices, close_file_line_indices):
                db2_file_path = lines[open_file_line_index].replace("open the file:", "").replace("Input ligand:", "").strip()
                db2_file_paths.append(db2_file_path)
                for data_row_line_index in range(open_file_line_index + 1, close_file_line_index):
                    data_row_line = lines[data_row_line_index]
                    data_row = data_row_line.strip().split()
                    if data_row[0].isdigit():
                        data_row = [db2_file_path] + data_row
                    else:
                        data_row = [db2_file_path] + [np.nan for _ in range(len(self.COLUMN_NAMES))]

                    # pad missing columns with NaN
                    if len(data_row) != len(df_column_names):
                        num_missing = len(df_column_names) - len(data_row)
                        data_row += [np.nan for _ in range(num_missing)]

                    #
                    data_row_dict = {column_name: data_row[i] for i, column_name in enumerate(df_column_names)}
                    data.append(data_row_dict)

            return pd.DataFrame.from_records(data)
