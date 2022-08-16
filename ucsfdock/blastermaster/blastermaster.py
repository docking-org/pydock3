import os
from dataclasses import fields
import logging

from ucsfdock.blastermaster.steps.receptor_most_occupied_residues_renaming import \
    ReceptorMostOccupiedResiduesRenamingStep
from ucsfdock.blastermaster.steps.receptor_protonation import ReceptorProtonationStep
from ucsfdock.blastermaster.steps.charged_receptor_deprotonation import ChargedReceptorDeprotonationStep
from ucsfdock.blastermaster.steps.ligand_hetatm_renaming import LigandHetatmRenamingStep
from ucsfdock.blastermaster.steps.binding_site_residues import BindingSiteResiduesSelectionStep
from ucsfdock.blastermaster.steps.molecular_surface import MolecularSurfaceGenerationStep
from ucsfdock.blastermaster.steps.binding_site_spheres import BindingSiteSpheresGenerationStep
from ucsfdock.blastermaster.steps.thin_spheres import ThinSpheresGenerationStep
from ucsfdock.blastermaster.steps.close_spheres import CloseSpheresGenerationStep
from ucsfdock.blastermaster.steps.pdb_to_sph import LigandPDBToSpheresConversionStep
from ucsfdock.blastermaster.steps.low_dielectric_spheres import LowDielectricSpheresSelectionStep
from ucsfdock.blastermaster.steps.sph_to_pdb import SpheresToPDBConversionStep
from ucsfdock.blastermaster.steps.matching_spheres import MatchingSpheresGenerationStep
from ucsfdock.blastermaster.steps.box import BoxGenerationStep
from ucsfdock.blastermaster.steps.receptor_transformation_for_electrostatics import \
    ReceptorTransformationForElectrostatics
from ucsfdock.blastermaster.steps.electrostatics import ElectrostaticsGridGenerationStepNoThinSpheres, \
    ElectrostaticsGridGenerationStepYesThinSpheres
from ucsfdock.blastermaster.steps.vdw import VDWScoringGridGenerationStep
from ucsfdock.blastermaster.steps.receptor_transformation_for_ligand_desolvation import \
    ReceptorTransformationForLigandDesolvationNoThinSpheres, ReceptorTransformationForLigandDesolvationYesThinSpheres
from ucsfdock.blastermaster.steps.ligand_desolvation import HydrogenAtomLigandDesolvationScoringGridGenerationStep, \
    HeavyAtomLigandDesolvationScoringGridGenerationStep
from ucsfdock.blastermaster.config import BlastermasterParametersConfiguration
from ucsfdock.util import get_dataclass_as_dict, get_logger_for_script
from ucsfdock.files import (
    Dir,
    File,
    IndockFile,
    INDOCK_FILE_NAME,
)
from ucsfdock.blastermaster.util import WorkingDir, BlasterFiles, BlasterFileNames
from ucsfdock.blastermaster import __file__ as BLASTERMASTER_INIT_FILE_PATH
from ucsfdock.blastermaster.defaults import __file__ as DEFAULTS_INIT_FILE_PATH


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#
BLASTER_TARGETS_DAG_PICKLE_FILE_NAME = "blaster_targets_dag.pickle"


def get_blaster_steps(blaster_files, param_dict, working_dir):
    #
    steps = []

    #
    def _get_step_dir(dir_name):
        return Dir(path=os.path.join(working_dir.path, dir_name))

    #
    steps.append(
        ReceptorMostOccupiedResiduesRenamingStep(
            step_dir=_get_step_dir("receptor_most_occupied_residues_renaming"),
            receptor_infile=blaster_files.receptor_file,
            receptor_most_occupied_residues_renamed_outfile=blaster_files.receptor_most_occupied_residues_renamed_file,
        )
    )

    #
    steps.append(
        ReceptorProtonationStep(
            step_dir=_get_step_dir("receptor_protonation"),
            receptor_infile=blaster_files.receptor_most_occupied_residues_renamed_file,
            add_h_dict_infile=blaster_files.add_h_dict_file,
            charged_receptor_outfile=blaster_files.charged_receptor_file,
            reduce_options_parameter=param_dict['receptor_protonation.reduce_options'],
        )
    )

    #
    if param_dict['covalent.use']:
        steps.append(
            ChargedReceptorDeprotonationStep(
                step_dir=_get_step_dir("charged_receptor_deprotonation"),
                charged_receptor_infile=blaster_files.charged_receptor_file,
                charged_receptor_deprotonated_outfile=blaster_files.charged_receptor_deprotonated_file,
                covalent_residue_num_parameter=param_dict['covalent.residue_num'],
                covalent_residue_name_parameter=param_dict['covalent.residue_name'],
                covalent_residue_atoms_parameter=param_dict['covalent.residue_atoms'],
            )
        )
        charged_receptor_file = blaster_files.charged_receptor_deprotonated_file
    else:
        charged_receptor_file = blaster_files.charged_receptor_file

    #
    steps.append(
        LigandHetatmRenamingStep(
            step_dir=_get_step_dir("ligand_hetatm_renaming"),
            ligand_infile=blaster_files.ligand_file,
            ligand_hetatm_renamed_outfile=blaster_files.ligand_hetatm_renamed_file,
        )
    )

    #
    steps.append(
        BindingSiteResiduesSelectionStep(
            step_dir=_get_step_dir("binding_site_residues_selection"),
            receptor_infile=blaster_files.receptor_file,
            ligand_infile=blaster_files.ligand_hetatm_renamed_file,
            filt_parameters_infile=blaster_files.binding_site_residues_parameters_file,
            binding_site_residues_outfile=blaster_files.binding_site_residues_file,
        )
    )

    #
    steps.append(
        MolecularSurfaceGenerationStep(
            step_dir=_get_step_dir("molecular_surface_generation"),
            charged_receptor_infile=charged_receptor_file,
            binding_site_residues_infile=blaster_files.binding_site_residues_file,
            radii_infile=blaster_files.molecular_surface_radii_file,
            molecular_surface_outfile=blaster_files.molecular_surface_file,
        )
    )

    #
    steps.append(
        BindingSiteSpheresGenerationStep(
            step_dir=_get_step_dir("binding_site_spheres_generation"),
            molecular_surface_infile=blaster_files.molecular_surface_file,
            spheres_outfile=blaster_files.all_spheres_file,
        )
    )

    #
    steps.append(
        LigandPDBToSpheresConversionStep(
            step_dir=_get_step_dir("ligand_pdb_to_spheres_conversion"),
            pdb_infile=blaster_files.ligand_hetatm_renamed_file,
            sph_outfile=blaster_files.ligand_matching_spheres_file,
        )
    )

    #
    steps.append(
        MatchingSpheresGenerationStep(
            step_dir=_get_step_dir("matching_spheres_generation"),
            charged_receptor_infile=charged_receptor_file,
            ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
            all_spheres_infile=blaster_files.all_spheres_file,
            matching_spheres_outfile=blaster_files.matching_spheres_file,
            covalent_use_parameter=param_dict['covalent.use'],
            covalent_residue_name_parameter=param_dict['covalent.residue_name'],
            covalent_residue_num_parameter=param_dict['covalent.residue_num'],
        )
    )

    #
    if param_dict['thin_spheres_elec.use']:
        #
        steps.append(
            MolecularSurfaceGenerationStep(
                step_dir=_get_step_dir("molecular_surface_generation_thin_spheres_elec"),
                charged_receptor_infile=charged_receptor_file,
                binding_site_residues_infile=blaster_files.binding_site_residues_file,
                radii_infile=blaster_files.molecular_surface_radii_file,
                molecular_surface_outfile=blaster_files.thin_spheres_elec_molecular_surface_file,
            )
        )

        #
        steps.append(
            ThinSpheresGenerationStep(
                step_dir=_get_step_dir("thin_spheres_generation_elec"),
                molecular_surface_infile=blaster_files.thin_spheres_elec_molecular_surface_file,
                thin_spheres_outfile=blaster_files.thin_spheres_elec_file,
                distance_to_surface_parameter=param_dict['thin_spheres_elec.distance_to_surface'],
                penetration_parameter=param_dict['thin_spheres_elec.penetration'],
            )
        )

        #
        steps.append(
            CloseSpheresGenerationStep(
                step_dir=_get_step_dir("close_spheres_generation_elec"),
                ligand_infile=blaster_files.ligand_hetatm_renamed_file,
                thin_spheres_infile=blaster_files.thin_spheres_elec_file,
                close_spheres_outfile=blaster_files.close_spheres_elec_file,
                distance_to_surface_parameter=param_dict['thin_spheres_elec.distance_to_surface'],
                penetration_parameter=param_dict['thin_spheres_elec.penetration'],
                distance_to_ligand_parameter=param_dict['thin_spheres_elec.distance_to_ligand'],
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                step_dir=_get_step_dir("close_spheres_elec_to_pdb_conversion"),
                sph_infile=blaster_files.close_spheres_elec_file,
                pdb_outfile=blaster_files.close_spheres_elec_pdb_file,
            )
        )

        #
        spheres_pdb_file_for_electrostatics = blaster_files.close_spheres_elec_pdb_file

    else:
        #
        steps.append(
            LowDielectricSpheresSelectionStep(
                step_dir=_get_step_dir("low_dielectric_spheres_selection"),
                charged_receptor_infile=charged_receptor_file,
                ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
                all_spheres_infile=blaster_files.all_spheres_file,
                low_dielectric_spheres_outfile=blaster_files.low_dielectric_spheres_file,
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                step_dir=_get_step_dir("low_dielectric_spheres_to_pdb_conversion"),
                sph_infile=blaster_files.low_dielectric_spheres_file,
                pdb_outfile=blaster_files.low_dielectric_spheres_pdb_file,
            )
        )

        #
        spheres_pdb_file_for_electrostatics = blaster_files.low_dielectric_spheres_pdb_file

    #
    if param_dict['thin_spheres_desolv.use']:
        #
        steps.append(
            MolecularSurfaceGenerationStep(
                step_dir=_get_step_dir("molecular_surface_generation_thin_spheres_desolv"),
                charged_receptor_infile=charged_receptor_file,
                binding_site_residues_infile=blaster_files.binding_site_residues_file,
                radii_infile=blaster_files.molecular_surface_radii_file,
                molecular_surface_outfile=blaster_files.thin_spheres_desolv_molecular_surface_file,
            )
        )

        #
        steps.append(
            ThinSpheresGenerationStep(
                step_dir=_get_step_dir("thin_spheres_generation_desolv"),
                molecular_surface_infile=blaster_files.thin_spheres_desolv_molecular_surface_file,
                thin_spheres_outfile=blaster_files.thin_spheres_desolv_file,
                distance_to_surface_parameter=param_dict['thin_spheres_desolv.distance_to_surface'],
                penetration_parameter=param_dict['thin_spheres_desolv.penetration'],
            )
        )

        #
        steps.append(
            CloseSpheresGenerationStep(
                step_dir=_get_step_dir("close_spheres_generation_desolv"),
                ligand_infile=blaster_files.ligand_hetatm_renamed_file,
                thin_spheres_infile=blaster_files.thin_spheres_desolv_file,
                close_spheres_outfile=blaster_files.close_spheres_desolv_file,
                distance_to_surface_parameter=param_dict['thin_spheres_desolv.distance_to_surface'],
                penetration_parameter=param_dict['thin_spheres_desolv.penetration'],
                distance_to_ligand_parameter=param_dict['thin_spheres_desolv.distance_to_ligand'],
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                step_dir=_get_step_dir("close_spheres_desolv_to_pdb_conversion"),
                sph_infile=blaster_files.close_spheres_desolv_file,
                pdb_outfile=blaster_files.close_spheres_desolv_pdb_file,
            )
        )

    #
    steps.append(
        BoxGenerationStep(
            step_dir=_get_step_dir("box_generation"),
            charged_receptor_infile=charged_receptor_file,
            ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
            box_outfile=blaster_files.box_file,
        )
    )

    #
    steps.append(
        ReceptorTransformationForElectrostatics(
            step_dir=_get_step_dir("receptor_transformation_for_electrostatics"),
            charged_receptor_infile=charged_receptor_file,
            spheres_pdb_infile=spheres_pdb_file_for_electrostatics,
            receptor_low_dielectric_pdb_outfile=blaster_files.receptor_low_dielectric_pdb_file,
        )
    )

    #
    if param_dict['thin_spheres_elec.use']:
        steps.append(
            ElectrostaticsGridGenerationStepYesThinSpheres(
                step_dir=_get_step_dir("electrostatics_grid_generation"),
                receptor_low_dielectric_pdb_infile=blaster_files.receptor_low_dielectric_pdb_file,
                charge_infile=blaster_files.electrostatics_charge_file,
                radius_infile=blaster_files.electrostatics_radius_file,
                delphi_infile=blaster_files.electrostatics_delphi_file,
                box_infile=blaster_files.box_file,
                electrostatics_phi_outfile=blaster_files.electrostatics_phi_file,
                electrostatics_pdb_outfile=blaster_files.electrostatics_pdb_file,
                electrostatics_trim_phi_outfile=blaster_files.electrostatics_trim_phi_file,
                electrostatics_phi_size_outfile=blaster_files.electrostatics_phi_size_file,
                thin_spheres_elec_distance_to_ligand_parameter=param_dict['thin_spheres_elec.distance_to_surface'],
                thin_spheres_elec_penetration_parameter=param_dict['thin_spheres_elec.penetration'],
            )
        )
    else:
        steps.append(
            ElectrostaticsGridGenerationStepNoThinSpheres(
                step_dir=_get_step_dir("electrostatics_grid_generation"),
                receptor_low_dielectric_pdb_infile=blaster_files.receptor_low_dielectric_pdb_file,
                charge_infile=blaster_files.electrostatics_charge_file,
                radius_infile=blaster_files.electrostatics_radius_file,
                delphi_infile=blaster_files.electrostatics_delphi_file,
                box_infile=blaster_files.box_file,
                electrostatics_phi_outfile=blaster_files.electrostatics_phi_file,
                electrostatics_pdb_outfile=blaster_files.electrostatics_pdb_file,
                electrostatics_trim_phi_outfile=blaster_files.electrostatics_trim_phi_file,
                electrostatics_phi_size_outfile=blaster_files.electrostatics_phi_size_file,
            )
        )

    #
    steps.append(
        VDWScoringGridGenerationStep(
            step_dir=_get_step_dir("vdw_scoring_grid_generation"),
            charged_receptor_infile=charged_receptor_file,
            vdw_parameters_infile=blaster_files.vdw_parameters_file,
            protein_table_infile=blaster_files.vdw_protein_table_file,
            box_infile=blaster_files.box_file,
            vdw_outfile=blaster_files.vdw_file,
            bump_map_outfile=blaster_files.vdw_bump_map_file,
        )
    )

    #
    if param_dict['thin_spheres_desolv.use']:
        steps.append(
            ReceptorTransformationForLigandDesolvationYesThinSpheres(
                step_dir=_get_step_dir("receptor_transformation_for_ligand_desolvation"),
                charged_receptor_pdb_infile=charged_receptor_file,
                close_spheres_desolv_pdb_infile=blaster_files.close_spheres_desolv_pdb_file,
                charged_receptor_desolv_pdb_outfile=blaster_files.charged_receptor_desolv_pdb_file,
            )
        )
    else:
        steps.append(
            ReceptorTransformationForLigandDesolvationNoThinSpheres(
                step_dir=_get_step_dir("receptor_transformation_for_ligand_desolvation"),
                charged_receptor_pdb_infile=charged_receptor_file,
                charged_receptor_desolv_pdb_outfile=blaster_files.charged_receptor_desolv_pdb_file,
            )
        )

    #
    steps.append(
        HydrogenAtomLigandDesolvationScoringGridGenerationStep(
            step_dir=_get_step_dir("ligand_desolvation_hydrogen"),
            box_infile=blaster_files.box_file,
            receptor_pdb_infile=blaster_files.charged_receptor_desolv_pdb_file,
            ligand_desolvation_outfile=blaster_files.ligand_desolvation_hydrogen_file,
            thin_spheres_desolv_use_parameter=param_dict['thin_spheres_desolv.use'],
            thin_spheres_desolv_distance_to_surface_parameter=param_dict['thin_spheres_desolv.distance_to_surface'],
            thin_spheres_desolv_penetration_parameter=param_dict['thin_spheres_desolv.penetration'],
            other_radius_parameter=param_dict['ligand_desolvation.other_radius'],
        )
    )

    #
    steps.append(
        HeavyAtomLigandDesolvationScoringGridGenerationStep(
            step_dir=_get_step_dir("ligand_desolvation_heavy"),
            box_infile=blaster_files.box_file,
            receptor_pdb_infile=blaster_files.charged_receptor_desolv_pdb_file,
            ligand_desolvation_outfile=blaster_files.ligand_desolvation_heavy_file,
            thin_spheres_desolv_use_parameter=param_dict['thin_spheres_desolv.use'],
            thin_spheres_desolv_distance_to_surface_parameter=param_dict['thin_spheres_desolv.distance_to_surface'],
            thin_spheres_desolv_penetration_parameter=param_dict['thin_spheres_desolv.penetration'],
            other_radius_parameter=param_dict['ligand_desolvation.other_radius'],
        )
    )

    return tuple(steps)


class Blastermaster(object):
    JOB_DIR_NAME = "blastermaster_job"
    CONFIG_FILE_NAME = "blastermaster_config.yaml"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.dirname(BLASTERMASTER_INIT_FILE_PATH), "default_blastermaster_config.yaml")
    WORKING_DIR_NAME = "working"
    DOCK_FILES_DIR_NAME = "dockfiles"
    DEFAULT_FILES_DIR_PATH = os.path.dirname(DEFAULTS_INIT_FILE_PATH)

    def __init__(self):
        pass
                
    def configure(self, job_dir_path=JOB_DIR_NAME, overwrite=False):
        # create job dir
        job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create working dir & copy in blaster files
        blaster_file_names = list(get_dataclass_as_dict(BlasterFileNames()).values())
        backup_blaster_file_paths = [os.path.join(self.DEFAULT_FILES_DIR_PATH, blaster_file_name) for blaster_file_name in blaster_file_names]
        blaster_file_names_in_cwd = [f for f in blaster_file_names if os.path.isfile(f)]
        files_to_copy_str = '\n\t'.join(blaster_file_names_in_cwd)
        if blaster_file_names_in_cwd:
            logger.info(f"Copying the following files from current directory into job working directory:\n\t{files_to_copy_str}")
        working_dir = WorkingDir(path=os.path.join(job_dir.path, self.WORKING_DIR_NAME), create=True, reset=False, files_to_copy_in=blaster_file_names_in_cwd, backup_files_to_copy_in=backup_blaster_file_paths)

        # create dock files dir
        dock_files_dir = Dir(path=os.path.join(job_dir.path, self.DOCK_FILES_DIR_NAME), create=True, reset=False)

        # write fresh config file from default file
        save_path = os.path.join(job_dir.path, self.CONFIG_FILE_NAME)
        BlastermasterParametersConfiguration.write_config_file(save_path, self.DEFAULT_CONFIG_FILE_PATH, overwrite=overwrite)

    def run(self, job_dir_path=".", config_file_path=None, use_graph_state=True, write_graph_image=False):
        # validate args
        if config_file_path is None:
            config_file_path = os.path.join(job_dir_path, self.CONFIG_FILE_NAME)
        try:
            File.validate_file_exists(config_file_path)
        except FileNotFoundError:
            logger.error("Config file not found. Are you in the job directory?")
            return

        # load directories
        job_dir = Dir(path=job_dir_path, create=True, reset=False)
        working_dir = WorkingDir(path=os.path.join(job_dir.path, self.WORKING_DIR_NAME), create=True, reset=False)
        dock_files_dir = Dir(path=os.path.join(job_dir.path, self.DOCK_FILES_DIR_NAME), create=True, reset=True)  # reset dock files dir in case re-running

        #
        blaster_files = BlasterFiles(working_dir=working_dir)
        indock_file = IndockFile(path=os.path.join(dock_files_dir.path, INDOCK_FILE_NAME))

        # load config file
        logger.info("Loading config file...")
        config = BlastermasterParametersConfiguration(config_file_path)
        logger.info("done.")

        # get params as str
        config_params_str = '\n'.join(
            [f"{param_name}: {param.value}" for param_name, param in config.param_dict.items()])
        logger.info(f"Parameters:\n{config_params_str}")

        # get blaster steps
        steps = get_blaster_steps(blaster_files, config.param_dict, working_dir)

        # reset step dirs
        for step in steps:
            if not step.is_done:
                step.step_dir.delete()

        # run steps
        logger.info("Running blaster steps in sequence.")
        for step in steps:
            step.run()

        # validate that all blaster files exist now
        # TODO

        # copy dock files to output directory
        logger.info("Copying dock files to dock files directory...")
        for dock_file_field in fields(blaster_files.dock_files):
            dock_file = getattr(blaster_files.dock_files, dock_file_field.name)
            File.copy_file(dock_file.path, os.path.join(dock_files_dir.path, dock_file.name))
        logger.info("done.")

        # make INDOCK file, using phi_size
        logger.info("Making indock file...")
        indock_file.write(blaster_files.dock_files, config.param_dict, dock_files_dir.name)
        logger.info("done.")


# TODO: flex blastermaster
'''
    def parse_res_list_list(in_list):
        """reads a structure in "55+56,67,89" and returns [[55,56],[67],[89]] """
        return [[int(id) for id in a_list.split("+")] for a_list in in_list.split(",")]

    def reform_res_list(res_list_list):
        """reads a list like [[55, 56], [67], [89]] and returns [55,56,67,89].
        yes, it removes the subgroupings."""
        return [residue for sub_list in res_list_list for residue in sub_list]

    def string_res_list(res_list):
        """reads a list like [55, 56] and returns "55+56" """
        return "+".join([str(residue).strip() for residue in res_list])

    def string_specific_char(combinations):
        """reads a set of combinations in
        [(number, occupancy, [residues], conformation),(...),(...)]
        and returns 1) "1.2.4.7" for instance, from the number and
        2) ["55A","56A","67B","89D"] combining residue and conformation"""

        # TODO: each combination should really be a namedtuple
        out_code = ".".join(
            [number for number, occupancy, residues, conformation in combinations]
        )
        out_res_list = [
            f"{residue}{conformation}"
            for number, occupancy, residues, conformation in combinations
            for residue in residues
        ]

        return out_code, out_res_list

    def make_flexible_files(self, groups, verbose=False):
        for combination in util.all_combinations(groups):
            out_code, out_res_list = self.string_specific_char(combination)

            if verbose:
                print(f"writing receptor combination file: {out_code}.pdb")

            pdb.specific_alts(self.blaster_files.receptor_file.path, out_res_list, f"{out_code}.pdb")

    def make_flexible_electrostatic_files(
        self,
        electrostatics_grid,
        groups,
        phi_size,
        scale_center,
        verbose=False,
    ):
        elec_dir = Dir("electrostatics_debugging_working", create=True, reset=False)  # TODO: put absolute path
        elec_dir_2 = Dir("electrostatics_debugging_working2", create=True, reset=False)  # TODO: put absolute path
        elec_dir_3 = Dir("elec_debug_working3", create=True, reset=False)  # for old-style elec. construct & compare # TODO: put absolute path
        elec_dir_out = Dir("electrostatics_debugging", create=True, reset=False)  # TODO: put absolute path

        # first do the old method of doing things
        movable_residues = list(
            set(
                [
                    res_list
                    for one_set_of_groups in groups
                    for group_num, occupancy, res_list, alt_char in one_set_of_groups
                ]
            )
        )

        # use atom charges in file, don't recalculate them,
        # neutralizes residues as desired
        scale_center_no_recharge = scale_center + ["input_atm_file=true"]
        for one_set_of_groups in groups:
            for group_num, occupancy, res_list, alt_char in one_set_of_groups:
                print(
                    group_num,
                    occupancy,
                    res_list,
                    alt_char,
                    qnifft_group_num_to_outfile_name_dict[group_num],
                )
                group_dir = os.path.join(elec_dir_3, str(group_num))
                os.mkdir(group_dir)
                pre_pdb = os.path.join(group_dir, "pre-neutral.pdb")
                post_pdb = os.path.join(group_dir, "post-neutral.pdb")
                shutil.copy(qnifft_group_num_to_outfile_name_dict[group_num], pre_pdb)
                pdb_data = pdb.PDBData(pre_pdb)
                if (
                    res_list == "invariant"
                ):  # special case, neutralize all but the reslist
                    new_pdb = pdb_data.clear_factors_residues(
                        movable_residues, matching=False
                    )
                else:  # neutralize the res_list
                    new_pdb = pdb_data.clear_factors_residues(res_list)
                new_pdb.write(post_pdb)
                self.programs.electrostatics_program.run(
                    self.log_files.electrostatics_log_file,
                    self.blaster_files.electrostatics_phi_file,
                    electrostatics_grid,
                    self.blaster_files.electrostatics_charge_file,
                    self.blaster_files.electrostatics_radius_file,
                    self.blaster_files.electrostatics_pdb_file,
                    "post-neutral.pdb",
                    group_dir,
                    False,
                    scale_center_no_recharge,
                )
                _, __ = phi.trim(
                    BlasterFile(os.path.join(group_dir, self.blaster_files.electrostatics_phi_file.name)),
                    self.blaster_files.box_file,
                    os.path.join(group_dir, self.blaster_files.electrostatics_trim_phi_file.name),
                )
        elec_dir.copy_in_file(
            self.blaster_files.low_dielectric_spheres_pdb_file.path,
            self.blaster_files.low_dielectric_spheres_pdb_file.name
        )
        for combination in util.all_combinations(groups):
            out_code, out_res_list = self.string_specific_char(combination)
            if verbose:
                print(f"writing receptor combination file: {out_code}.pdb")
            pdb.specific_alts(
                self.blaster_files.receptor_file.path,
                out_res_list,
                os.path.join(elec_dir, f"{out_code}.temp.pdb"),
            )
            pdb.delete_alt_chars(
                os.path.join(elec_dir, f"{out_code}.temp.pdb"),
                os.path.join(elec_dir, f"{out_code}.postalt.temp.pdb"),
            )
            self.programs.add_hydrogens_program.run(
                working_dir=elec_dir,
                receptor_infile=f"{out_code}.postalt.temp.pdb",
                add_h_dict_infile=self.blaster_files.add_h_dict_file,
                charged_receptor_outfile=f"{out_code}.pdb",
            )
            shutil.copy(
                os.path.join(elec_dir, f"{out_code}.pdb"),
                os.path.join(elec_dir_out, f"{out_code}.pdb"),
            )
            cat(
                infile_1=File(os.path.join(elec_dir, f"{out_code}.pdb")),
                infile_2=File(os.path.join(elec_dir, self.blaster_files.low_dielectric_spheres_pdb_file)),
                outfile=File(os.path.join(elec_dir, f"{out_code}.sph.pdb")),
            )

            self.programs.electrostatics_program.run(
                f"{out_code}{electrostatics_log_file_path}",  # TODO
                f"{out_code}{electrostatics_outfile_path}",  # TODO
                electrostatics_grid,
                electrostatics_charge_file_path,
                electrostatics_radius_file_path,
                out_code + electrostatics_pdb_outfile_path,
                f"{out_code}.sph.pdb",
                elec_dir,
                False,
                scale_center,
            )
            phi_size_check, phi_center_check = phi.trim(
                os.path.join(
                    elec_dir, out_code + self.blaster_files.electrostatics_phi_file.name
                ),
                self.blaster_files.box_file.path,
                os.path.join(elec_dir, out_code + self.blaster_files.electrostatics_trim_phi_file.name),
            )
            shutil.copy(
                os.path.join(elec_dir, out_code + self.blaster_files.electrostatics_trim_phi_file.name),
                os.path.join(elec_dir_out, out_code + self.blaster_files.electrostatics_trim_phi_file.name),
            )
            out_code_list = out_code.split(".")
            # compute sums of phimaps to emulate the DOCK score of any atom
            working_phi_sum = os.path.join(elec_dir_2, f"{out_code}.sum.phi")
            working_phi_sum_temp = os.path.join(elec_dir_2, f"{out_code}.temp.sum.phi")
            phi.add(
                os.path.join(
                    outfile_dir_path,
                    str(out_code_list[0]),
                    electrostatics_trim_outfile_path,
                ),
                os.path.join(
                    outfile_dir_path,
                    str(out_code_list[1]),
                    electrostatics_trim_outfile_path,
                ),
                working_phi_sum,
                phi_size,
            )
            for out_code_one in out_code_list[2:]:
                shutil.copy(working_phi_sum, working_phi_sum_temp)
                phi.add(
                    working_phi_sum_temp,
                    os.path.join(
                        outfile_dir_path,
                        str(out_code_one),
                        electrostatics_trim_outfile_path,
                    ),
                    working_phi_sum,
                    phi_size,
                )
            shutil.copy(
                working_phi_sum, os.path.join(elec_dir_out, f"{out_code}.sum.phi")
            )
            phi.subtract(
                os.path.join(elec_dir_out, out_code + self.blaster_files.electrostatics_trim_phi_file.name),
                os.path.join(elec_dir_out, f"{out_code}.sum.phi"),
                os.path.join(elec_dir_out, f"{out_code}.diff.phi"),
                phi_size,
            )

            # same thing but for old style grids
            working_phi_sum = os.path.join(elec_dir_3, out_code + ".sum.phi")
            working_phi_sum_temp = os.path.join(elec_dir_3, out_code + ".temp.sum.phi")
            phi.add(
                os.path.join(
                    elec_dir_3, str(out_code_list[0]), self.blaster_files.electrostatics_trim_phi_file.name
                ),
                os.path.join(
                    elec_dir_3, str(out_code_list[1]), self.blaster_files.electrostatics_trim_phi_file.name
                ),
                working_phi_sum,
                phi_size,
            )
            for out_code_one in out_code_list[2:]:
                shutil.copy(working_phi_sum, working_phi_sum_temp)
                phi.add(
                    working_phi_sum_temp,
                    os.path.join(
                        elec_dir_3, str(out_code_one), self.blaster_files.electrostatics_trim_phi_file.name
                    ),
                    working_phi_sum,
                    phi_size,
                )
            shutil.copy(
                working_phi_sum, os.path.join(elec_dir_out, f"{out_code}.old.sum.phi")
            )
            phi.subtract(
                os.path.join(
                    elec_dir_out, f"{out_code}{self.blaster_files.electrostatics_trim_phi_file.name}"
                ),
                os.path.join(elec_dir_out, f"{out_code}.old.sum.phi"),
                os.path.join(elec_dir_out, f"{out_code}.old.diff.phi"),
                phi_size,
            )

    def make_flexible_readme(self, groups, flexible_readme, verbose=False):
        """writes a file explaining the filenames/codes that are used to represent
        all possible flexible receptor combinations"""
        if verbose:
            print(
                "writing flexible receptor explanation file:",
                flexible_readme,
            )
        with open(flexible_readme, "w") as out_file:
            for combination in util.all_combinations(groups):
                code = ".".join(
                    [str(flexible_group[0]) for flexible_group in combination]
                )
                out_file.write(f"{code} ")
                out_file.write(f"{combination[0][2]} {combination[0][3]} ")
                for flexible_group in combination[1:]:
                    out_file.write(f"{self.string_res_list(flexible_group[2])} ")
                    out_file.write(f"{str(flexible_group[3])} ")
                out_file.write("\n")

    def make_part_readme(self, groups, covalent_readme, verbose=False):
        """writes a file explaining the filenames/codes that are used to represent
        each part. include energy"""
        if verbose:
            print("writing part explanation file:", covalent_readme)
        with open(covalent_readme, "w") as out_file:
            for one_set in groups:
                for grids_dir_count, occupancy, one_res_list, alt_char in one_set:
                    out_file.write(str(grids_dir_count))
                    if one_res_list == "invariant":
                        out_file.write(" invariant")
                    else:
                        out_file.write(" " + self.string_res_list(one_res_list))
                    out_file.write(" " + alt_char)
                    out_file.write(
                        " "
                        + str(util.occupancy_to_energy(occupancy, flexible_penalty_m))
                    )
                    out_file.write(" " + str(occupancy))
                    out_file.write("\n")

    def run_flex(self):
        """rewrite of MakeDOCK/DOCKblaster into non-makefile form. this file controls
        everything about the whole process. see default input files in
        blast_defaultfiles. based on previous versions by Irwin, Mysinger, Lorber,
        Wei, Kirschner, and Huang.
        this version makes flexible receptor grid files, which was never done
        by a previous version."""

        # TODO: this whole function needs its filepath variable manipulation to be reworked

        outfile_paths = []

        # just beginning, lots more files get added to this by the end
        output_0_files = []  # implicitly 0 files, used in electrostatics when there
        # really is no change between the invariant pose and this pose in the
        # most occupied position

        # reset working dir
        util.reset_dir(self.working_dir_path)

        # copy input files to working dir
        receptor_file_path = util.copy_file(
            self.blaster_files_config.receptor_infile_name, self.working_dir_path
        )
        ligand_file_path = util.copy_file(ligand_infile_path, self.working_dir_path)
        util.copy_file(
            os.path.join(util.DEFAULTS_DIR_PATH, vdw_parameters),
            self.working_dir_path,
        )
        dock_files_dirs = []  # empty now, will have directories
        most_occupied_pdb = "most.occupied.pdb"
        if verbose:
            print(
                "copying input files into working directory, fixing columns",
                receptor_infile_path,
                ligand_infile_path,
                self.working_dir_path,
            )
        decoded_residue_lists = self.parse_res_list_list(flexible_residues)
        all_flex_residues = self.reform_res_list(decoded_residue_lists)
        if add_hydrogens_first:
            unprot_receptor = "noprot." + receptor_file_path
            temp_receptor = "tempprot." + receptor_file_path
            pdb.move_columns(
                receptor_file_path,
                os.path.join(self.working_dir.path, unprot_receptor),
            )
            add_hydrogens(
                working_dir_path=self.working_dir_path,
                add_h_program_file_path=add_h_program_file_path,
                receptor_infile_name=unprot_receptor,
                add_h_dict_infile_name=add_h_dict_file_name,
                charged_receptor_outfile_name=temp_receptor,
                add_h_options="-OH -HIS -ALLALT -ROTNH3 -Keep",
            )

            pdb.del_hydrogens(
                os.path.join(self.working_dir.path, temp_receptor),
                os.path.join(self.working_dir.path, receptor_file_path),
                all_flex_residues,
            )
        else:  # just copy unprotonated to receptor file
            pdb.move_columns(
                receptor_file_path,
                receptor_file_path,  # TODO: make regular file_path
            )
        pdb.move_columns(
            ligand_file_path,
            ligand_file_path,  # TODO: make regular file_path
        )
        add_hydrogens(
            working_dir_path=working_dir_path,
            add_h_program_file_path=add_h_program_file_path,
            receptor_infile_name=receptor_file_name,
            add_h_dict_infile_name=add_h_dict_file_name,
            charged_receptor_outfile_name=add_h_outfile_name,
        )
        # first deviation, call most occupied to get a receptor with just one position
        pdb.most_occupied(
            os.path.join(self.working_dir.path, receptor_file_path),
            os.path.join(self.working_dir.path, most_occupied_pdb),
        )
        most_occupied_del_alt_char_pdb = "most.occupied.nochar.pdb"
        pdb.delete_alt_chars(
            os.path.join(self.working_dir.path, most_occupied_pdb),
            os.path.join(self.working_dir.path, most_occupied_del_alt_char_pdb),
        )
        add_hydrogens(
            working_dir_path=self.working_dir_path,
            add_h_program_file_path=add_h_program_file_path,
            receptor_infile_name=most_occupied_del_alt_char_pdb,
            add_h_dict_infile_name=add_h_dict_file_name,
            charged_receptor_outfile_name=add_h_outfile_name,
        )
        # lots of tasks in this subfunction
        qnifft_pdb_out_names = {}  # used for debugging electrostatics
        (
            phi_size,
            matching_spheres_file_path,
            electrostatics_trim_phi_file_path,
        ) = self.run_single_receptor(
            util.BlasterFiles.HYDROGEN_ADDITION_OUTFILE_NAME,
            os.path.join(util.DEFAULTS_DIR_PATH, "radii"),
        )

        outfile_paths.append(matching_spheres_file_path)
        outfile_paths.append(electrostatics_trim_phi_file_path)

        qnifft_pdb_out_names[1] = os.path.join(
            self.working_dir_path, electrostatics_pdb_output
        )
        # figure out the scale and center for qnifft run, will use for all qnifft runs
        scale_center = get_scale_and_center(self.log_files.electrostatics_log_file)
        non_flex_dir = "1"  # non-flexible part always = 1
        groups = [[(1, 1.0, "invariant", "-")]]
        # list of lists, where each sublist is
        # the set of possibilities for a choice between many alternates.
        # the number is the occupancy, the text is the flexible residues
        non_flex_working_dir = Dir(os.path.join(self.working_dir.path, non_flex_dir))
        dock_files_dirs.append(non_flex_dir)
        non_flex_working_dir.copy_in_file(self.blaster_files.electrostatics_trim_phi_file.path)
        outfile_paths.append(os.path.join(non_flex_dir, self.blaster_files.electrostatics_trim_phi_file.name))
        # make extra files for flexible receptor stuff.
        # this stuff is for the invariant receptor that has all residues with alter.
        # conformations deleted entirely. for vdw and ligand ligand_desolvation only.
        pdb.delete_alts(
            self.blaster_files.charged_receptor_file.path,
            os.path.join(
                self.working_dir.path,
                non_flex_dir,
                self.blaster_files.charged_receptor_file.name,
            ),
            all_flex_residues,
        )
        self.blaster_files.box_file.copy_to(
            os.path.join(
                self.working_dir.path, non_flex_dir, self.blaster_files.box_file.name
            )
        )
        shutil.copy(
            os.path.join(
                self.working_dir.path,
                self.blaster_files.thin_spheres_desolv_file + ".close.pdb",
            ),
            os.path.join(
                self.working_dir.path,
                non_flex_dir,
                self.blaster_files.thin_spheres_desolv_file + ".close.pdb",
            ),
        )
        self.run_single_receptor_just_vdw_ligand_desolv(
            os.path.join(self.working_dir.path, non_flex_dir),
            copy_to_dir=os.path.join(self.dock_files_dir.path, non_flex_dir),
        )
        # these files need saved
        outfile_paths.append(os.path.join(non_flex_dir, self.blaster_files.vdw_bump_map_file))
        outfile_paths.append(os.path.join(non_flex_dir, self.blaster_files.vdw_file))
        outfile_paths.append(
            os.path.join(
                non_flex_dir,
                self.blaster_files.ligand_desolvation_hydrogen_file.name,
            )
        )
        outfile_paths.append(
            os.path.join(
                non_flex_dir,
                self.blaster_files.ligand_desolvation_heavy_file.name,
            )
        )
        # now proceed to the flexible parts
        alt_prefix = "alts.full"
        grids_dir_count = 2
        pdb.move_columns(
            receptor_file_path,
            os.path.join(self.working_dir.path, "temp." + receptor_file_path),
        )
        orig_pdb_data = pdb.PDBData(
            os.path.join(self.working_dir.path, "temp." + receptor_file_path)
        )
        for one_res_list in decoded_residue_lists:  # do each set separately
            alt_names, alt_chars = pdb.make_alts(
                os.path.join(self.working_dir.path, receptor_file_path),
                os.path.join(self.working_dir.path, alt_prefix),
                [one_res_list],
            )
            this_group = []
            for i, alt_name in enumerate(alt_names):
                # print "XXX", alt_count, alt_name, alt_chars[alt_count], one_res_list,
                is_this_most_occupied = orig_pdb_data.is_most_occupied_residue_chain(
                    one_res_list[0], alt_chars[i]
                )
                # print is_this_most_occupied
                alt_name_before_no_char = alt_name + ".bnc.pdb"
                alt_name_before_no_char_pre_h = alt_name + ".bnc.prepolarh.pdb"
                shutil.move(alt_name, alt_name_before_no_char)
                pdb.delete_alt_chars(
                    alt_name_before_no_char, alt_name_before_no_char_pre_h
                )
                add_hydrogens(
                    working_dir_path=self.working_dir_path,
                    add_h_program_file_path=add_h_program_file_path,
                    receptor_infile_name=os.path.basename(
                        alt_name_before_no_char_pre_h
                    ),
                    add_h_dict_infile_name=add_h_dict_file_name,
                    charged_receptor_outfile_name=os.path.basename(
                        os.path.basename(alt_name)
                    ),
                )
                pdb_data = pdb.PDBData(alt_name)  # to get occupancy data
                occupancy = pdb_data.get_occupancy_residue(one_res_list[0])
                alt_char = alt_chars[i]
                # make directory to store everything in
                grids_dir = os.path.join(self.working_dir.path, str(grids_dir_count))
                this_group.append((grids_dir_count, occupancy, one_res_list, alt_char))
                os.mkdir(grids_dir)
                dock_files_dirs.append(str(grids_dir_count))
                # temporary names used -- could be user controllable at some point
                temp_name = "rec.justpart.pdb"
                sub_phi_name = "subtracted.phi"
                # copy files
                shutil.copy(alt_name, os.path.join(grids_dir, temp_name))
                shutil.copy(
                    os.path.join(
                        self.working_dir_path,
                        sph_to_pdb_low_dielectric_pdb_output,
                    ),
                    grids_dir,
                )
                shutil.copy(
                    os.path.join(
                        self.working_dir_path,
                        thin_spheres_sph_output_elec + ".close.pdb",
                    ),
                    grids_dir,
                )

                # do phimap making stuff
                if thin_spheres_use:
                    cat(
                        infile_path_1=os.path.join(grids_dir, temp_name),
                        infile_path_2=os.path.join(
                            grids_dir, thin_spheres_sph_output_elec + ".close.pdb"
                        ),
                        outfile_path=os.path.join(
                            grids_dir,
                            receptor_with_charges_receptor_low_dielectric_output,
                        ),
                    )
                else:
                    cat(
                        infile_path_1=os.path.join(grids_dir, temp_name),
                        infile_path_2=os.path.join(
                            grids_dir, sph_to_pdb_low_dielectric_pdb_output
                        ),
                        outfile_path=os.path.join(
                            grids_dir,
                            receptor_with_charges_receptor_low_dielectric_output,
                        ),
                    )

                qnifft_pdb_out_names[grids_dir_count] = os.path.join(
                    grids_dir, electrostatics_pdb_output
                )
                if (
                    not is_this_most_occupied
                ) or debug_electrostatics:  # TODO: Ask someone if this logic is correct
                    # if this is most occupied, electrostatics are 0
                    electrostatics(
                        util.BlasterPrograms.ELECTROSTATICS_PROGRAM_PATH,
                        self.log_files_config.electrostatics_log_name,
                        util.BlasterFiles.ELECTROSTATICS_OUTFILE_NAME,
                        electrostatics.grid,
                        electrostatics_charge_file,
                        electrostatics_radius_file,
                        electrostatics_pdb_output,
                        receptor_with_charges_receptor_low_dielectric_output,
                        grids_dir,
                        False,
                        scale_center,
                    )
                    phi.subtract(
                        os.path.join(
                            grids_dir, util.BlasterFiles.ELECTROSTATICS_OUTFILE_NAME
                        ),
                        os.path.join(
                            self.working_dir_path,
                            util.BlasterFiles.ELECTROSTATICS_OUTFILE_NAME,
                        ),
                        os.path.join(grids_dir, sub_phi_name),
                    )
                    phi_size_check, phi_center_check = phi.trim(
                        os.path.join(grids_dir, sub_phi_name),
                        os.path.join(
                            self.working_dir_path, util.BlasterFiles.BOX_OUTFILE_NAME
                        ),
                        os.path.join(grids_dir, electrostatics_trim_output),
                    )
                    if phi_size != phi_size_check:  # should be the same
                        print("something wrong, trimmed phimaps are not the same size.")
                        print(phi_size, phi_size_check)
                        sys.exit(1)
                    # add to files necessary for docking
                    outfile_paths.append(
                        os.path.join(str(grids_dir_count), electrostatics_trim_output)
                    )
                else:  # we don't need to do electrostatics, implicitly 0 everywhere
                    output_0_files.append(
                        (str(grids_dir_count), "electrostatics")
                    )  # tuple of
                    # grid# and then what grid is 0 (electrostatics only one for now,
                    # more likely added later like recdes and hydrophobic effect
                # now do vdw & ligand desolv
                # copy files
                grids_del_dir = grids_dir + ".del"
                os.mkdir(grids_del_dir)
                alt_del_name = alt_name + ".del.pdb"
                pdb.del_all_but(alt_name, alt_del_name, one_res_list)
                shutil.copy(
                    alt_del_name,
                    os.path.join(
                        grids_del_dir, util.BlasterFiles.HYDROGEN_ADDITION_OUTFILE_NAME
                    ),
                )
                shutil.copy(
                    alt_del_name,
                    os.path.join(grids_del_dir, self.blaster_files.charged_receptor_desolv_pdb_file),
                )
                shutil.copy(
                    os.path.join(
                        self.working_dir_path, util.BlasterFiles.BOX_OUTFILE_NAME
                    ),
                    os.path.join(grids_del_dir, util.BlasterFiles.BOX_OUTFILE_NAME),
                )
                # make an empty sph file
                print(
                    "touch "
                    + os.path.join(
                        grids_del_dir, self.blaster_files.thin_spheres_desolv_file + ".close.pdb"
                    )
                )
                os.system(
                    "touch "
                    + os.path.join(
                        grids_del_dir, self.blaster_files.thin_spheres_desolv_file + ".close.pdb"
                    )
                )
                # run, but make sure scripts put the files in the right place if sge
                self.run_single_receptor_just_vdw_ligand_desolv(
                    grids_del_dir,
                    copy_to_dir=os.path.join(
                        self.dock_files_dir_path, str(grids_dir_count)
                    ),
                )
                # these files need copied and saved
                shutil.copy(
                    os.path.join(grids_del_dir, vdw_prefix + ".bmp"),
                    os.path.join(grids_dir, vdw_prefix + ".bmp"),
                )
                outfile_paths.append(
                    os.path.join(str(grids_dir_count), vdw_prefix + ".bmp")
                )
                shutil.copy(
                    os.path.join(grids_del_dir, vdw_prefix + ".vdw"),
                    os.path.join(grids_dir, vdw_prefix + ".vdw"),
                )
                outfile_paths.append(
                    os.path.join(str(grids_dir_count), vdw_prefix + ".vdw")
                )
                try:
                    shutil.copy(
                        os.path.join(
                            grids_del_dir,
                            ligand_desolvation_prefix
                            + ligand_desolvation_hydrogen_name,
                        ),
                        os.path.join(
                            grids_dir,
                            ligand_desolvation_prefix
                            + ligand_desolvation_hydrogen_name,
                        ),
                    )
                except IOError:
                    pass  # error handled later
                outfile_paths.append(
                    os.path.join(
                        str(grids_dir_count),
                        ligand_desolvation_prefix + ligand_desolvation_hydrogen_name,
                    )
                )
                try:
                    shutil.copy(
                        os.path.join(
                            grids_del_dir,
                            ligand_desolvation_prefix + ligand_desolvation_heavy_name,
                        ),
                        os.path.join(
                            grids_dir,
                            ligand_desolvation_prefix + ligand_desolvation_heavy_name,
                        ),
                    )
                except IOError:
                    pass  # error handled later
                outfile_paths.append(
                    os.path.join(
                        str(grids_dir_count),
                        ligand_desolvation_prefix + ligand_desolvation_heavy_name,
                    )
                )
                # advance counter of grid numbered directories
                grids_dir_count += 1
            groups.append(this_group)
        # want to copy all files into output directory
        if not os.path.exists(self.dock_files_dir_path):
            os.mkdir(self.dock_files_dir_path)  # make working directory if it isn't there
            for sub_dir in dock_files_dirs:
                os.mkdir(os.path.join(self.dock_files_dir_path, sub_dir))
        for outfile_path in outfile_paths:
            print(f"copying {outfile_path} into {self.dock_files_dir_path}")
            sub_dir = os.path.dirname(outfile_path)
            try:
                shutil.copy(
                    os.path.join(self.working_dir.path, outfile_path),
                    os.path.join(self.dock_files_dir_path, sub_dir),
                )
            except IOError:  # file missing
                print(
                    "docking file MISSING:",
                    os.path.join(self.working_dir.path, outfile_path),
                )
        # make INDOCK file
        make_indock(phi_size, "..", groups, output_0_files)
        # make files that represent each possible recombination of receptors
        self.make_flexible_files(receptor_file_path, groups)
        # make readme files
        self.make_flexible_readme(groups)
        self.make_part_readme(groups)
        # make flexible electrostatic files for debugging electrostatic approximations
        if debug_electrostatics:
            self.make_flexible_electrostatic_files(
                self.working_dir_path,
                receptor_file_path,
                groups,
                phi_size,
                scale_center,
                qnifft_pdb_out_names,
            )
'''

# TODO: blaster state file
'''
        # load blaster targets graph from previous run if it exists
        if self.use_graph_state:
            self.blaster_targets_dag_pickle_file_path = os.path.join(self.job_dir.path, BLASTER_TARGETS_DAG_PICKLE_FILE_NAME)
            if File.file_exists(self.blaster_targets_dag_pickle_file_path):
                with open(self.blaster_targets_dag_pickle_file_path, 'rb') as f:
                    blaster_targets_dag_previous_run = pickle.load(f)

                # compare blaster target graphs and validate
                if not nx.algorithms.isomorphism.is_isomorphic(self.blaster_targets_dag, blaster_targets_dag_previous_run):
                    raise Exception(
                        "Blaster targets graph for this run and blaster targets graph for previous run are not isomorphic!")
                for node_name in self.blaster_targets_dag.nodes:
                    try:
                        blaster_target_current_run = self.blaster_targets_dag.nodes[node_name]["blaster_target"]
                        blaster_target_previous_run = blaster_targets_dag_previous_run.nodes[node_name]["blaster_target"]
                        assert type(blaster_target_current_run) == type(blaster_target_previous_run)
                    except:
                        raise Exception(
                            "Blaster target graphs for current run and previous run are of different types!")

                #
                nodes_edited = []
                nodes_deleted = []
                for node_name in self.blaster_targets_dag.nodes:
                    blaster_target_current_run = self.blaster_targets_dag.nodes[node_name]["blaster_target"]
                    blaster_target_previous_run = blaster_targets_dag_previous_run.nodes[node_name]["blaster_target"]
                    if blaster_target_current_run.exists:
                        if blaster_target_previous_run.datetime_marked_complete_in_most_recent_job is not None:
                            if blaster_target_current_run.datetime_marked_complete_in_most_recent_job > blaster_target_previous_run.datetime_marked_complete_in_most_recent_job:
                                nodes_edited.append(node_name)
                    else:
                        if blaster_target_previous_run.datetime_marked_complete_in_most_recent_job is not None:
                            nodes_deleted.append(node_name)
                logger.info(f"Blaster targets edited since previous job: {nodes_edited}")
                logger.info(f"Blaster targets deleted since previous job: {nodes_deleted}")

                #
                nodes_whose_non_edited_descendants_have_already_been_deleted = []
                nodes_edited_set = set(nodes_edited)

                def _delete_non_edited_descendants(node_name):
                    if node_name in nodes_whose_non_edited_descendants_have_already_been_deleted:
                        return

                    if node_name not in nodes_edited_set:
                        if node_name in nodes_dependent_on_config_parameter_to_be_optimized:
                            self.sub_jobs_dir.delete()
                        self.blaster_targets_dag.nodes[node_name]["blaster_target"].delete()
                        nodes_whose_non_edited_descendants_have_already_been_deleted.append(node_name)

                    successors = list(self.blaster_targets_dag.successors(node_name))
                    if successors:
                        for successor in successors:
                            _delete_non_edited_descendants(successor)

                #
                for node_name in nodes_edited + nodes_deleted:
                    _delete_non_edited_descendants(node_name)

            #
            if self.write_graph_image:
                nx.draw_networkx(self.blaster_targets_dag, arrowsize=20, node_size=500,
                                 pos=nx.spring_layout(self.blaster_targets_dag,
                                                      k=7 / math.sqrt(self.blaster_targets_dag.order())))
                fig = plt.gcf()
                fig.set_size_inches(20, 20)
                plt.savefig(os.path.join(self.job_dir.path, "blaster_targets_graph.png"))
'''
