import os
import logging
import re
import collections
from copy import deepcopy
from datetime import datetime
from functools import wraps
from dataclasses import make_dataclass

from pydock3.util import validate_variable_type, system_call
from pydock3.config import Parameter
from pydock3.files import File, Dir, ProgramFile, LogFile
from pydock3.blastermaster.programs import __file__ as PROGRAMS_INIT_FILE_PATH
from pydock3.blastermaster.defaults import __file__ as DEFAULTS_INIT_FILE_PATH


#
PROGRAMS_DIR_PATH = os.path.dirname(PROGRAMS_INIT_FILE_PATH)
DEFAULT_FILES_DIR_PATH = os.path.dirname(DEFAULTS_INIT_FILE_PATH)

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT = {
    'add_h_dict_file': "reduce_wwPDB_het_dict.txt",
    'binding_site_residues_parameters_file': "filt.params",
    'molecular_surface_radii_file': "radii",
    'electrostatics_charge_file': "amb.crg.oxt",
    'electrostatics_radius_file': "vdw.siz",
    'electrostatics_delphi_file': "delphi.def",
    'vdw_parameters_file': "vdw.parms.amb.mindock",
    'vdw_protein_table_file': "prot.table.ambcrg.ambH",
    'residue_code_to_polar_h_yaml_file': "residue_code_polar_h.yaml",
    'receptor_file': "rec.pdb",
    'ligand_file': "xtal-lig.pdb",
    'receptor_most_occupied_residues_renamed_file': "rec.most_occ_renamed.pdb",
    'ligand_hetatm_renamed_file': "xtal-lig.hetatm_renamed.pdb",
    'charged_receptor_file': "rec.crg.pdb",
    'charged_receptor_deprotonated_file': "rec.crg.deprotonated.pdb",
    'binding_site_residues_file': "rec.site",
    'molecular_surface_file': "rec.ms",
    'all_spheres_file': "all_spheres.sph",
    'electrostatics_phi_file': "qnifft.electrostatics.phi",
    'electrostatics_pdb_file': "qnifft.atm",
    'electrostatics_trim_phi_file': "trim.electrostatics.phi",
    'electrostatics_phi_size_file': "phi.size",
    'thin_spheres_elec_file': "thin_spheres_elec.sph",
    'thin_spheres_desolv_file': "thin_spheres_desolv.sph",
    'thin_spheres_elec_molecular_surface_file': "rec.ts_elec.ms",
    'thin_spheres_desolv_molecular_surface_file': "rec.ts_desolv.ms",
    'close_spheres_elec_file': f"thin_spheres_elec.sph.close",
    'close_spheres_desolv_file': f"thin_spheres_desolv.sph.close",
    'close_spheres_elec_pdb_file': f"thin_spheres_elec.sph.close.pdb",
    'close_spheres_desolv_pdb_file': f"thin_spheres_desolv.close.pdb",
    'matching_spheres_file': "matching_spheres.sph",
    'box_file': "box",
    'ligand_matching_spheres_file': "xtal-lig.match.sph",
    'low_dielectric_spheres_file': "lowdielectric.sph",
    'low_dielectric_spheres_pdb_file': f"lowdielectric.sph.pdb",
    'receptor_low_dielectric_pdb_file': "receptor.crg.lowdielectric.pdb",
    'charged_receptor_desolv_pdb_file': "rec.crg.lds.pdb",
    'vdw_file': "vdw.vdw",
    'vdw_bump_map_file': "vdw.bmp",
    'ligand_desolvation_heavy_file': "ligand.desolv.heavy",
    'ligand_desolvation_hydrogen_file': "ligand.desolv.hydrogen",
}
PROPER_BLASTER_FILE_NAME_TO_BLASTER_FILE_IDENTIFIER_DICT = {value: key for key, value in BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.items()}
DOCK_FILE_IDENTIFIERS = [
    "electrostatics_phi_size_file",
    "electrostatics_trim_phi_file",
    "ligand_desolvation_heavy_file",
    "ligand_desolvation_hydrogen_file",
    "matching_spheres_file",
    "vdw_bump_map_file",
    "vdw_file",
    "vdw_parameters_file",
]
DOCK_FILE_IDENTIFIER_TO_PROPER_DOCK_FILE_NAME_DICT = {dock_file_identifier: BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT[dock_file_identifier] for dock_file_identifier in DOCK_FILE_IDENTIFIERS}


#
class ProgramFilePaths:
    REDUCE_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "reduce/reduce")
    DMS_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "dms/bin/dms")
    FILT_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "filt/bin/filt")
    SPHGEN_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "sphgen/bin/sphgen")
    PDBTOSPH_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "pdbtosph/bin/pdbtosph"
    )
    MAKESPHERES1_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "makespheres1/makespheres1.cli.pl"
    )
    DOSHOWSPH_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "showsphere/doshowsph.csh"
    )
    MAKESPHERES3_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "makespheres3/makespheres3.cli.pl"
    )
    QNIFFT_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "qnifft/bin/qnifft22_193_pgf_32"
    )
    MAKEBOX_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "makebox/makebox.smallokay.pl"
    )
    CHEMGRID_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "chemgrid/bin/chemgrid"
    )
    SOLVMAP_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "solvmap/bin/solvmap")


class BlasterFile(File):
    """#TODO"""

    def __init__(self, path, identifier, src_file_path=None):
        super().__init__(path=path)

        #
        self.identifier = identifier

        #
        if src_file_path is not None:
            if self.exists:
                self.src_file_path = None
                logger.debug(
                    f"Source file provided for {self.name} but file already exists in working dir: '{self.path}'. Using existing version instead of source file {src_file_path}."
                )
            else:
                self.src_file_path = src_file_path
                self.copy_from(self.src_file_path, overwrite=True)
        else:
            self.src_file_path = None

        #
        if self.exists:
            self.datetime_marked_complete_in_most_recent_job = (
                self.get_datetime_file_was_last_modified(self.path)
            )
        else:
            self.datetime_marked_complete_in_most_recent_job = None

    @property
    def src_file_path(self):
        return self._src_file_path

    @src_file_path.setter
    def src_file_path(self, src_file_path):
        if src_file_path is not None:
            File.validate_file_exists(src_file_path)
        self._src_file_path = src_file_path

    @property
    def datetime_marked_complete_in_most_recent_job(self):
        return self._datetime_marked_complete_in_most_recent_job

    @datetime_marked_complete_in_most_recent_job.setter
    def datetime_marked_complete_in_most_recent_job(self, dt):
        if self.exists:
            validate_variable_type(dt, allowed_instance_types=(datetime,))
        else:
            validate_variable_type(dt, allowed_instance_types=(type(None),))
        self._datetime_marked_complete_in_most_recent_job = dt

    def copy_from(self, src_file_path, overwrite=True):
        self.copy_file(
            src_file_path=src_file_path, dst_file_path=self.path, overwrite=overwrite
        )
        self.datetime_marked_complete_in_most_recent_job = (
            self.get_datetime_file_was_last_modified(self.path)
        )

    def delete(self):
        self.delete_file(self.path)
        self.src_file_path = None
        self.datetime_marked_complete_in_most_recent_job = None

    def __eq__(self, other):
        if type(other) == type(self):
            return (self.path == other.path) and (self.identifier == other.identifier)
        return False


class WorkingDir(Dir):
    """#TODO"""

    def __init__(
        self,
        path,
        create=False,
        reset=False,
        files_to_copy_in=None,
        new_file_names=None,
        backup_files_to_copy_in=None,
        new_backup_file_names=None,
    ):
        super().__init__(path, create=create, reset=reset)

        #
        if files_to_copy_in is None:
            files_to_copy_in = []
        if new_file_names is None:
            new_file_names = []
        if backup_files_to_copy_in is None:
            backup_files_to_copy_in = []
        if new_backup_file_names is None:
            new_backup_file_names = []

        #
        if len(files_to_copy_in) != len(new_file_names):
            raise Exception("# files to copy in must match # of new file names.")
        if len(backup_files_to_copy_in) != len(new_backup_file_names):
            raise Exception("# backup files to copy in must match # of new backup file names.")

        # copy in specified files if they exist, otherwise try to copy in backup files
        file_names_to_copy_in = [
            File.get_file_name_of_file(file_path) for file_path in files_to_copy_in
        ]
        for src_backup_file_path, dst_backup_file_name in zip(
            backup_files_to_copy_in, new_backup_file_names
        ):
            if (
                File.get_file_name_of_file(src_backup_file_path)
                not in file_names_to_copy_in
            ):
                if File.file_exists(src_backup_file_path):
                    self.copy_in_file(
                        src_backup_file_path, dst_file_name=dst_backup_file_name
                    )
        for src_file_path, dst_file_name in zip(files_to_copy_in, new_file_names):
            if File.file_exists(src_file_path):
                self.copy_in_file(src_file_path, dst_file_name=dst_file_name)


class BlasterFiles(object):
    """See: https://wiki.bkslab.org/index.php/Blastermaster_files"""

    def __init__(
        self,
        working_dir,
    ):
        #
        for blaster_file_identifier, proper_blaster_file_name in BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.items():
            blaster_file = BlasterFile(path=os.path.join(working_dir.path, proper_blaster_file_name), identifier=blaster_file_identifier)
            setattr(self, blaster_file_identifier, blaster_file)

    @property
    def dock_files(self):
        return DockFiles(**{dock_file_identifier: getattr(self, dock_file_identifier) for dock_file_identifier in DOCK_FILE_IDENTIFIERS})

    def get_attribute_name_of_blaster_file_with_file_name(self, file_name):
        attributes = [
            a
            for a in dir(self)
            if not a.startswith("__") and not callable(getattr(self, a))
        ]
        for a in attributes:
            if isinstance(getattr(self, a), (BlasterFile,)):
                if getattr(self, a).name == file_name:
                    return a

        raise Exception(
            f"Blaster file with file name '{file_name}' not found in attributes of BlasterFiles instance."
        )


DockFiles = make_dataclass("DockFiles", [(identifier, BlasterFile) for identifier in DOCK_FILE_IDENTIFIERS])


class BlasterStep(object):
    def __init__(self, working_dir, infile_tuples, outfile_tuples, parameter_tuples, program_file_path=None):
        #
        self.step_dir = self._get_step_dir(working_dir, outfile_tuples)

        #
        self._infiles = None
        self._outfiles = None
        self._parameters = None

        #
        self.infiles = self._process_infiles(*infile_tuples)
        self.outfiles = self._process_outfiles(*outfile_tuples)
        self.parameters = self._process_parameters(*parameter_tuples)

        #
        if program_file_path:
            self.program_file = ProgramFile(program_file_path)
        else:
            self.program_file = None

        #
        self.log_file = self._get_log_file()

    def __str__(self):
        return self.__class__.__name__

    def _get_step_dir(self, working_dir, outfile_tuples):
        class_name_snake_case = re.sub('(?<!^)(?=[A-Z])', '_', self.__str__()).lower()
        comma_separated_outfile_names = ','.join([x[0].name for x in outfile_tuples])
        dir_name = f"{class_name_snake_case}_outfiles={comma_separated_outfile_names}"
        return Dir(path=os.path.join(working_dir.path, dir_name))

    @property
    def infiles(self):
        return self._infiles

    @infiles.setter
    def infiles(self, infiles_named_tuple):
        validate_variable_type(infiles_named_tuple, allowed_instance_types=(tuple,))
        for obj in infiles_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(BlasterFile,))

        self._infiles = infiles_named_tuple

    def _process_infiles(self, *infile_tuples):
        step_infiles = []
        for i, (infile, arg_name, new_file_name) in enumerate(infile_tuples):
            #
            validate_variable_type(infile, allowed_instance_types=(BlasterFile,))
            File.validate_path(infile.path)
            validate_variable_type(arg_name, allowed_instance_types=(str,))
            validate_variable_type(new_file_name, allowed_instance_types=(str, type(None),))

            #
            step_infile = deepcopy(infile)

            #
            if new_file_name is not None:
                step_infile.path = os.path.join(
                    self.step_dir.path, new_file_name
                )
            else:
                step_infile.path = os.path.join(self.step_dir.path, infile.name)
            step_infile.original_file_in_working_dir = infile

            #
            step_infiles.append(step_infile)

        #
        Infiles = collections.namedtuple(
            "Infiles", " ".join([arg_name for infile, arg_name, new_file_name in infile_tuples])
        )

        #
        return Infiles(*step_infiles)

    @property
    def outfiles(self):
        return self._outfiles

    @outfiles.setter
    def outfiles(self, outfiles_named_tuple):
        validate_variable_type(outfiles_named_tuple, allowed_instance_types=(tuple,))
        for obj in outfiles_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(BlasterFile,))

        self._outfiles = outfiles_named_tuple

    def _process_outfiles(self, *outfile_tuples):
        step_outfiles = []
        for i, (outfile, arg_name, new_file_name) in enumerate(outfile_tuples):
            #
            validate_variable_type(outfile, allowed_instance_types=(BlasterFile,))
            File.validate_path(outfile.path)
            validate_variable_type(arg_name, allowed_instance_types=(str,))
            validate_variable_type(new_file_name, allowed_instance_types=(str, type(None),))

            #
            step_outfile = deepcopy(outfile)

            #
            if new_file_name is not None:
                step_outfile.path = os.path.join(
                    self.step_dir.path, new_file_name
                )
            else:
                step_outfile.path = os.path.join(self.step_dir.path, outfile.name)
            step_outfile.original_file_in_working_dir = outfile

            #
            step_outfiles.append(step_outfile)

        #
        Outfiles = collections.namedtuple(
            "Outfiles", " ".join([arg_name for outfile, arg_name, new_file_name in outfile_tuples])
        )

        #
        return Outfiles(*step_outfiles)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, parameters_named_tuple):
        validate_variable_type(parameters_named_tuple, allowed_instance_types=(tuple,))
        for obj in parameters_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(Parameter,))

        self._parameters = parameters_named_tuple

    def _process_parameters(self, *parameter_tuples):
        step_parameters = []
        for parameter, arg_name in parameter_tuples:
            #
            validate_variable_type(parameter, allowed_instance_types=(Parameter,))
            validate_variable_type(arg_name, allowed_instance_types=(str,))

            #
            step_parameter = deepcopy(parameter)

            #
            step_parameters.append(step_parameter)

        #
        Parameters = collections.namedtuple(
            "Parameters",
            " ".join([arg_name for parameter, arg_name in parameter_tuples]),
        )

        #
        return Parameters(*step_parameters)

    @property
    def step_dir(self):
        return self._step_dir

    @step_dir.setter
    def step_dir(self, step_dir):
        Dir.validate_obj_is_dir(step_dir)
        self._step_dir = step_dir

    @property
    def is_done(self):
        #
        if self.infiles is None or self.outfiles is None:
            logger.exception(
                f"Variables 'infiles' and 'outfiles' not defined for {self.__class__.__name__} instance."
            )
            raise

        #
        return all(
            [outfile.original_file_in_working_dir.exists for outfile in self.outfiles]
        )

    @staticmethod
    def handle_run_func(run_func):
        @wraps(run_func)
        def wrapper(self):
            if self.is_done:
                logger.info(f"Skipping {self.__class__.__name__} since is_done=True")
            else:
                logger.info(f"Running {self.__class__.__name__}")
                self._set_up_step_dir()
                run_func(self)
                self._export_outfiles()

        return wrapper

    def _set_up_step_dir(self):
        self.step_dir.create(reset=True)
        self._import_infiles()

    def run(self):
        raise NotImplementedError

    def _get_log_file(self):
        return LogFile(path=os.path.join(self.step_dir.path, "log"))

    def _import_infiles(self):
        """Copies input files to be used in program process into step dir"""
        for infile in self.infiles:
            infile.copy_from(infile.original_file_in_working_dir.path)

    def _export_outfiles(self):
        """Copies output files created in program process back out into working dir"""
        for outfile in self.outfiles:
            if not outfile.original_file_in_working_dir.exists:
                outfile.original_file_in_working_dir.copy_from(outfile.path)

    def run_command(self, command_str, timeout_seconds=None, env_vars_dict=None):
        result = system_call(
            command_str,
            cwd=self.step_dir.path,
            timeout_seconds=timeout_seconds,
            env_vars_dict=env_vars_dict,
        )
        with open(self.log_file.path, "a") as f:
            f.write(f"command:\n{command_str}\n")
            f.write("\n")
            f.write(f"stdout:\n{result.stdout}\n")
            f.write("\n")
            f.write(f"stderr:\n{result.stderr}\n")
            f.write("\n")
            f.write("-" * 20)
            f.write("\n\n")

    def log_parameters_file(self, file):
        logger.debug(
            f"{self.program_file.name} parameters file at {file.path} begins: "
        )
        with open(file.path, "r") as f:
            logger.debug("\n\t".join([line.strip() for line in f.readlines()]))
        logger.debug(f"\t{self.program_file.name} parameters file ends.")
