import os
import logging
import collections
from copy import deepcopy
from datetime import datetime
from functools import wraps
from dataclasses import dataclass

from pydock3.util import validate_variable_type, system_call
from pydock3.config import Parameter
from pydock3.files import File, Dir, LogFile
from pydock3.blastermaster.programs import __file__ as PROGRAMS_INIT_FILE_PATH

PROGRAMS_DIR_PATH = os.path.dirname(PROGRAMS_INIT_FILE_PATH)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


#
class ProgramFilePaths:
    REDUCE_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "reduce/reduce")
    DMS_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "dms/bin/dms")
    FILT_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "filt/bin/filt")
    SPHGEN_PROGRAM_FILE_PATH = os.path.join(PROGRAMS_DIR_PATH, "sphgen/bin/sphgen")
    THIN_SPHERES_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "thinspheres/thin_spheres.py"
    )
    CLOSE_SPH_PROGRAM_FILE_PATH = os.path.join(
        PROGRAMS_DIR_PATH, "thinspheres/close_sph.py"
    )
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

    def __init__(self, path, src_file_path=None):
        super().__init__(path=path)

        #
        self.src_file_path = None

        #

        if src_file_path is not None:
            if self.exists:
                logger.debug(f"Source file provided for {self.name} but file already exists in working dir: '{self.path}'. Using existing version instead of source file {src_file_path}.")
            else:
                self.src_file_path = src_file_path
                self.copy_from(self.src_file_path, overwrite=True)

        #
        if self.exists:
            self.datetime_marked_complete_in_most_recent_job = self.get_datetime_file_was_last_modified(self.path)
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
        self.copy_file(src_file_path=src_file_path, dst_file_path=self.path, overwrite=overwrite)
        self.datetime_marked_complete_in_most_recent_job = self.get_datetime_file_was_last_modified(self.path)

    def delete(self):
        self.delete_file(self.path)
        self.src_file_path = None
        self.datetime_marked_complete_in_most_recent_job = None


class WorkingDir(Dir):
    """#TODO"""

    def __init__(self, path, create=False, reset=False, files_to_copy_in=None, backup_files_to_copy_in=None):
        super().__init__(path, create=create, reset=reset)

        #
        if files_to_copy_in is None:
            files_to_copy_in = []
        if backup_files_to_copy_in is None:
            backup_files_to_copy_in = []

        # copy in specified files if they exist, otherwise try to copy in backup files
        file_names_to_copy_in = [File.get_file_name_of_file(file_path) for file_path in files_to_copy_in]
        for backup_file_path in backup_files_to_copy_in:
            if File.get_file_name_of_file(backup_file_path) not in file_names_to_copy_in:
                if File.file_exists(backup_file_path):
                    self.copy_in_file(backup_file_path)
        for file_path in files_to_copy_in:
            if File.file_exists(file_path):
                self.copy_in_file(file_path)


@dataclass
class BlasterFileNames(object):
    #
    add_h_dict_file_name: str = "reduce_wwPDB_het_dict.txt"
    binding_site_residues_parameters_file_name: str = "filt.params"
    molecular_surface_radii_file_name: str = "radii"
    electrostatics_charge_file_name: str = "amb.crg.oxt"
    electrostatics_radius_file_name: str = "vdw.siz"
    electrostatics_delphi_file_name: str = "delphi.def"
    vdw_parameters_file_name: str = "vdw.parms.amb.mindock"
    vdw_protein_table_file_name: str = "prot.table.ambcrg.ambH"

    #
    receptor_file_name: str = "rec.pdb"
    ligand_file_name: str = "xtal-lig.pdb"

    receptor_most_occupied_residues_renamed_file_name: str = "rec.most_occ_renamed.pdb"

    ligand_hetatm_renamed_file_name: str = "xtal-lig.hetatm_renamed.pdb"

    charged_receptor_file_name: str = "rec.crg.pdb"

    charged_receptor_deprotonated_file_name: str = "rec.crg.deprotonated.pdb"

    binding_site_residues_file_name: str = "rec.site"

    molecular_surface_file_name: str = "rec.ms"

    all_spheres_file_name: str = "all_spheres.sph"

    electrostatics_phi_file_name: str = "qnifft.electrostatics.phi"
    electrostatics_pdb_file_name: str = "qnifft.atm"
    electrostatics_trim_phi_file_name: str = "trim.electrostatics.phi"
    electrostatics_phi_size_file_name: str = "phi.size"

    thin_spheres_elec_file_name: str = "thin_spheres_elec.sph"
    thin_spheres_desolv_file_name: str = "thin_spheres_desolv.sph"

    thin_spheres_elec_molecular_surface_file_name: str = "rec.ts_elec.ms"
    thin_spheres_desolv_molecular_surface_file_name: str = "rec.ts_desolv.ms"

    close_spheres_elec_file_name: str = f"{thin_spheres_elec_file_name}.close"
    close_spheres_desolv_file_name: str = f"{thin_spheres_desolv_file_name}.close"

    close_spheres_elec_pdb_file_name: str = f"{thin_spheres_elec_file_name}.close.pdb"
    close_spheres_desolv_pdb_file_name: str = f"{thin_spheres_desolv_file_name}.close.pdb"

    matching_spheres_file_name: str = "matching_spheres.sph"

    box_file_name: str = "box"

    ligand_matching_spheres_file_name: str = "xtal-lig.match.sph"

    low_dielectric_spheres_file_name: str = "lowdielectric.sph"
    low_dielectric_spheres_pdb_file_name: str = f"{low_dielectric_spheres_file_name}.pdb"

    receptor_low_dielectric_pdb_file_name: str = "receptor.crg.lowdielectric.pdb"

    charged_receptor_desolv_pdb_file_name: str = "rec.crg.lds.pdb"

    vdw_file_name: str = "vdw.vdw"
    vdw_bump_map_file_name: str = "vdw.bmp"

    ligand_desolvation_heavy_file_name: str = "ligand.desolv.heavy"
    ligand_desolvation_hydrogen_file_name: str = "ligand.desolv.hydrogen"


class BlasterFiles(object):
    """See: https://wiki.bkslab.org/index.php/Blastermaster_files"""

    def __init__(
            self,
            working_dir,
    ):
        #
        blaster_file_names = BlasterFileNames()

        #
        self.add_h_dict_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.add_h_dict_file_name))
        self.binding_site_residues_parameters_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.binding_site_residues_parameters_file_name))
        self.molecular_surface_radii_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.molecular_surface_radii_file_name))
        self.electrostatics_charge_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_charge_file_name))
        self.electrostatics_radius_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_radius_file_name))
        self.electrostatics_delphi_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_delphi_file_name))
        self.vdw_parameters_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.vdw_parameters_file_name))
        self.vdw_protein_table_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.vdw_protein_table_file_name))
        self.receptor_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.receptor_file_name))
        self.receptor_most_occupied_residues_renamed_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.receptor_most_occupied_residues_renamed_file_name))
        self.ligand_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.ligand_file_name))
        self.ligand_hetatm_renamed_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.ligand_hetatm_renamed_file_name))
        self.charged_receptor_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.charged_receptor_file_name))
        self.charged_receptor_deprotonated_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.charged_receptor_deprotonated_file_name))
        self.binding_site_residues_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.binding_site_residues_file_name))
        self.molecular_surface_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.molecular_surface_file_name))
        self.all_spheres_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.all_spheres_file_name))
        self.electrostatics_phi_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_phi_file_name))
        self.electrostatics_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_pdb_file_name))
        self.electrostatics_trim_phi_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_trim_phi_file_name))
        self.electrostatics_phi_size_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.electrostatics_phi_size_file_name))
        self.thin_spheres_elec_molecular_surface_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.thin_spheres_elec_molecular_surface_file_name))
        self.thin_spheres_desolv_molecular_surface_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.thin_spheres_desolv_molecular_surface_file_name))
        self.thin_spheres_elec_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.thin_spheres_elec_file_name))
        self.thin_spheres_desolv_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.thin_spheres_desolv_file_name))
        self.close_spheres_elec_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.close_spheres_elec_file_name))
        self.close_spheres_desolv_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.close_spheres_desolv_file_name))
        self.close_spheres_elec_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.close_spheres_elec_pdb_file_name))
        self.close_spheres_desolv_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.close_spheres_desolv_pdb_file_name))
        self.matching_spheres_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.matching_spheres_file_name))
        self.box_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.box_file_name))
        self.ligand_matching_spheres_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.ligand_matching_spheres_file_name))
        self.low_dielectric_spheres_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.low_dielectric_spheres_file_name))
        self.low_dielectric_spheres_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.low_dielectric_spheres_pdb_file_name))
        self.receptor_low_dielectric_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.receptor_low_dielectric_pdb_file_name))
        self.charged_receptor_desolv_pdb_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.charged_receptor_desolv_pdb_file_name))
        self.vdw_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.vdw_file_name))
        self.vdw_bump_map_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.vdw_bump_map_file_name))
        self.ligand_desolvation_heavy_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.ligand_desolvation_heavy_file_name))
        self.ligand_desolvation_hydrogen_file = BlasterFile(path=os.path.join(working_dir.path, blaster_file_names.ligand_desolvation_hydrogen_file_name))

    @property
    def dock_files(self):
        return DockFiles(
            matching_spheres_file=self.matching_spheres_file,
            electrostatics_trim_phi_file=self.electrostatics_trim_phi_file,
            vdw_file=self.vdw_file,
            vdw_bump_map_file=self.vdw_bump_map_file,
            vdw_parameters_file=self.vdw_parameters_file,
            ligand_desolvation_heavy_file=self.ligand_desolvation_heavy_file,
            ligand_desolvation_hydrogen_file=self.ligand_desolvation_hydrogen_file,
            electrostatics_phi_size_file=self.electrostatics_phi_size_file,
        )

    def get_attribute_name_of_blaster_file_with_file_name(self, file_name):
        attributes = [a for a in dir(self) if not a.startswith('__') and not callable(getattr(self, a))]
        for a in attributes:
            if isinstance(getattr(self, a), (BlasterFile,)):
                if getattr(self, a).name == file_name:
                    return a

        logger.exception(f"Blaster file with file name '{file_name}' not found in attributes of BlasterFiles instance.")
        raise


@dataclass
class DockFiles:
    matching_spheres_file: BlasterFile
    electrostatics_trim_phi_file: BlasterFile
    vdw_file: BlasterFile
    vdw_bump_map_file: BlasterFile
    vdw_parameters_file: BlasterFile
    ligand_desolvation_heavy_file: BlasterFile
    ligand_desolvation_hydrogen_file: BlasterFile
    electrostatics_phi_size_file: BlasterFile


class BlasterStep(object):
    def __init__(self, step_dir):
        #
        self.step_dir = step_dir

        #
        self.program_file = None
        self._infiles = None
        self._outfiles = None
        self._parameters = None

        #
        self.log_file = self._get_log_file()

    def __str__(self):
        return self.__class__.__name__

    @property
    def infiles(self):
        return self._infiles

    @infiles.setter
    def infiles(self, infiles_named_tuple):
        validate_variable_type(infiles_named_tuple, allowed_instance_types=(tuple,))
        for obj in infiles_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(BlasterFile,))

        self._infiles = infiles_named_tuple

    def process_infiles(self, *infile_tuples, new_file_names_tuple=None):
        step_infiles = []
        for i, (infile, arg_name) in enumerate(infile_tuples):
            #
            validate_variable_type(infile, allowed_instance_types=(BlasterFile,))
            File.validate_path(infile.path)

            #
            step_infile = deepcopy(infile)

            #
            if new_file_names_tuple is not None:
                step_infile.path = os.path.join(self.step_dir.path, new_file_names_tuple[i])
            else:
                step_infile.path = os.path.join(self.step_dir.path, infile.name)
            step_infile.original_file_in_working_dir = infile

            #
            step_infiles.append(step_infile)

        #
        Infiles = collections.namedtuple("Infiles", " ".join([arg_name for infile, arg_name in infile_tuples]))

        #
        self.infiles = Infiles(*step_infiles)

    @property
    def outfiles(self):
        return self._outfiles

    @outfiles.setter
    def outfiles(self, outfiles_named_tuple):
        validate_variable_type(outfiles_named_tuple, allowed_instance_types=(tuple,))
        for obj in outfiles_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(BlasterFile,))

        self._outfiles = outfiles_named_tuple

    def process_outfiles(self, *outfile_tuples, new_file_names_tuple=None):
        step_outfiles = []
        for i, (outfile, arg_name) in enumerate(outfile_tuples):
            #
            validate_variable_type(outfile, allowed_instance_types=(BlasterFile,))
            File.validate_path(outfile.path)

            #
            step_outfile = deepcopy(outfile)

            #
            if new_file_names_tuple is not None:
                step_outfile.path = os.path.join(self.step_dir.path, new_file_names_tuple[i])
            else:
                step_outfile.path = os.path.join(self.step_dir.path, outfile.name)
            step_outfile.original_file_in_working_dir = outfile

            #
            step_outfiles.append(step_outfile)

        #
        Outfiles = collections.namedtuple("Outfiles", " ".join([arg_name for outfile, arg_name in outfile_tuples]))

        #
        self.outfiles = Outfiles(*step_outfiles)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, parameters_named_tuple):
        validate_variable_type(parameters_named_tuple, allowed_instance_types=(tuple,))
        for obj in parameters_named_tuple:
            validate_variable_type(obj, allowed_instance_types=(Parameter,))

        self._parameters = parameters_named_tuple

    def process_parameters(self, *parameter_tuples):
        step_parameters = []
        for parameter, arg_name in parameter_tuples:
            #
            validate_variable_type(parameter, allowed_instance_types=(Parameter,))

            #
            step_parameter = deepcopy(parameter)

            #
            step_parameters.append(step_parameter)

        #
        Parameters = collections.namedtuple("Parameters", " ".join([arg_name for parameter, arg_name in parameter_tuples]))

        #
        self.parameters = Parameters(*step_parameters)

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
        return all([outfile.original_file_in_working_dir.exists for outfile in self.outfiles])

    @staticmethod
    def handle_run_func(run_func):
        @wraps(run_func)
        def wrapper(self):
            if self.is_done:
                logger.info(
                    f"Skipping {self.__class__.__name__} since is_done=True"
                )
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
        result = system_call(command_str, cwd=self.step_dir.path, timeout_seconds=timeout_seconds, env_vars_dict=env_vars_dict)
        with open(self.log_file.path, "a") as f:
            f.write(f"command:\n{command_str}\n")
            f.write("\n")
            f.write(f"stdout:\n{result.stdout}\n")
            f.write("\n")
            f.write(f"stderr:\n{result.stderr}\n")
            f.write("\n")
            f.write("-"*20)
            f.write("\n\n")

    def log_parameters_file(self, file):
        logger.debug(
            f"{self.program_file.name} parameters file at {file.path} begins: "
        )
        with open(file.path, "r") as f:
            logger.debug("\n\t".join([line.strip() for line in f.readlines()]))
        logger.debug(f"\t{self.program_file.name} parameters file ends.")
