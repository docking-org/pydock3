import os
from dataclasses import fields, astuple
import logging

from pydock3.blastermaster.steps.receptor_most_occupied_residues_renaming import (
    ReceptorMostOccupiedResiduesRenamingStep,
)
from pydock3.blastermaster.steps.receptor_protonation import ReceptorProtonationStep
from pydock3.blastermaster.steps.charged_receptor_deprotonation import (
    ChargedReceptorDeprotonationStep,
)
from pydock3.blastermaster.steps.ligand_hetatm_renaming import LigandHetatmRenamingStep
from pydock3.blastermaster.steps.binding_site_residues import (
    BindingSiteResiduesSelectionStep,
)
from pydock3.blastermaster.steps.molecular_surface import MolecularSurfaceGenerationStep
from pydock3.blastermaster.steps.binding_site_spheres import (
    BindingSiteSpheresGenerationStep,
)
from pydock3.blastermaster.steps.thin_spheres import ThinSpheresGenerationStep
from pydock3.blastermaster.steps.close_spheres import CloseSpheresGenerationStep
from pydock3.blastermaster.steps.pdb_to_sph import LigandPDBToSpheresConversionStep
from pydock3.blastermaster.steps.low_dielectric_spheres import (
    LowDielectricSpheresSelectionStep,
)
from pydock3.blastermaster.steps.sph_to_pdb import SpheresToPDBConversionStep
from pydock3.blastermaster.steps.matching_spheres import MatchingSpheresGenerationStep
from pydock3.blastermaster.steps.box import BoxGenerationStep
from pydock3.blastermaster.steps.receptor_transformation_for_electrostatics import (
    ReceptorTransformationForElectrostatics,
)
from pydock3.blastermaster.steps.electrostatics import (
    ElectrostaticsGridGenerationStepNoThinSpheres,
    ElectrostaticsGridGenerationStepYesThinSpheres,
)
from pydock3.blastermaster.steps.vdw import VDWScoringGridGenerationStep
from pydock3.blastermaster.steps.receptor_transformation_for_ligand_desolvation import (
    ReceptorTransformationForLigandDesolvationNoThinSpheres,
    ReceptorTransformationForLigandDesolvationYesThinSpheres,
)
from pydock3.blastermaster.steps.ligand_desolvation import (
    HydrogenAtomLigandDesolvationScoringGridGenerationStep,
    HeavyAtomLigandDesolvationScoringGridGenerationStep,
)
from pydock3.blastermaster.steps.visualization import VisualizationStep
from pydock3.blastermaster.config import BlastermasterParametersConfiguration
from pydock3.util import Script, get_dataclass_as_dict, unpack_step_params
from pydock3.config import flatten_and_parameter_cast_param_dict
from pydock3.files import (
    Dir,
    File,
    IndockFile,
    INDOCK_FILE_NAME,
)
from pydock3.blastermaster.util import BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT, WorkingDir, BlasterFiles
from pydock3.blastermaster import __file__ as BLASTERMASTER_INIT_FILE_PATH
from pydock3.blastermaster.defaults import __file__ as DEFAULTS_INIT_FILE_PATH


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#
BLASTER_TARGETS_DAG_PICKLE_FILE_NAME = "blaster_targets_dag.pickle"
BINARY_DOCK_FILE_IDENTIFIERS = {"vdw_file", "electrostatics_trim_phi_file"}


def get_blaster_steps(blaster_files, flat_param_dict, working_dir):
    #
    steps = []

    #
    steps.append(
        ReceptorMostOccupiedResiduesRenamingStep(
            working_dir=working_dir,
            receptor_infile=blaster_files.receptor_file,
            receptor_most_occupied_residues_renamed_outfile=blaster_files.receptor_most_occupied_residues_renamed_file,
        )
    )

    #
    steps.append(
        ReceptorProtonationStep(
            working_dir=working_dir,
            receptor_infile=blaster_files.receptor_most_occupied_residues_renamed_file,
            add_h_dict_infile=blaster_files.add_h_dict_file,
            residue_code_polar_h_yaml_infile=blaster_files.residue_code_to_polar_h_yaml_file,
            charged_receptor_outfile=blaster_files.charged_receptor_file,
            reduce_options_parameter=flat_param_dict[
                "receptor_protonation.reduce_options"
            ],
        )
    )

    #
    if flat_param_dict["covalent.use"]:
        steps.append(
            ChargedReceptorDeprotonationStep(
                working_dir=working_dir,
                charged_receptor_infile=blaster_files.charged_receptor_file,
                charged_receptor_deprotonated_outfile=blaster_files.charged_receptor_deprotonated_file,
                covalent_residue_num_parameter=flat_param_dict["covalent.residue_num"],
                covalent_residue_name_parameter=flat_param_dict["covalent.residue_name"],
                covalent_residue_atoms_parameter=flat_param_dict["covalent.residue_atoms"],
            )
        )
        blaster_files.charged_receptor_file = (blaster_files.charged_receptor_deprotonated_file)

    #
    steps.append(
        LigandHetatmRenamingStep(
            working_dir=working_dir,
            ligand_infile=blaster_files.ligand_file,
            ligand_hetatm_renamed_outfile=blaster_files.ligand_hetatm_renamed_file,
        )
    )

    #
    steps.append(
        BindingSiteResiduesSelectionStep(
            working_dir=working_dir,
            receptor_infile=blaster_files.charged_receptor_file,
            ligand_infile=blaster_files.ligand_hetatm_renamed_file,
            filt_parameters_infile=blaster_files.binding_site_residues_parameters_file,
            binding_site_residues_outfile=blaster_files.binding_site_residues_file,
        )
    )

    #
    steps.append(
        MolecularSurfaceGenerationStep(
            working_dir=working_dir,
            charged_receptor_infile=blaster_files.charged_receptor_file,
            binding_site_residues_infile=blaster_files.binding_site_residues_file,
            radii_infile=blaster_files.molecular_surface_radii_file,
            molecular_surface_outfile=blaster_files.molecular_surface_file,
            **(
                {'molecular_surface_density_parameter': flat_param_dict['molecular_surface_density']}
                if 'molecular_surface_density' in flat_param_dict
                else {}
            )
        )
    )

    #
    steps.append(
        BindingSiteSpheresGenerationStep(
            working_dir=working_dir,
            molecular_surface_infile=blaster_files.molecular_surface_file,
            spheres_outfile=blaster_files.all_spheres_file,
        )
    )

    #
    steps.append(
        LigandPDBToSpheresConversionStep(
            working_dir=working_dir,
            pdb_infile=blaster_files.ligand_hetatm_renamed_file,
            sph_outfile=blaster_files.ligand_matching_spheres_file,
        )
    )

    #
    steps.append(
        MatchingSpheresGenerationStep(
            working_dir=working_dir,
            charged_receptor_infile=blaster_files.charged_receptor_file,
            ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
            all_spheres_infile=blaster_files.all_spheres_file,
            matching_spheres_outfile=blaster_files.matching_spheres_file,
            covalent_use_parameter=flat_param_dict["covalent.use"],
            covalent_residue_name_parameter=flat_param_dict["covalent.residue_name"],
            covalent_residue_num_parameter=flat_param_dict["covalent.residue_num"],
            **unpack_step_params(flat_param_dict, "matching_spheres_generation")
        )
    )

    #
    if flat_param_dict["thin_spheres_elec.use"]:
        #
        steps.append(
            MolecularSurfaceGenerationStep(
                working_dir=working_dir,
                charged_receptor_infile=blaster_files.charged_receptor_file,
                binding_site_residues_infile=blaster_files.binding_site_residues_file,
                radii_infile=blaster_files.molecular_surface_radii_file,
                molecular_surface_outfile=blaster_files.thin_spheres_elec_molecular_surface_file,
                molecular_surface_density_parameter=flat_param_dict["thin_spheres_elec.molecular_surface_density"],
            )
        )

        #
        steps.append(
            ThinSpheresGenerationStep(
                working_dir=working_dir,
                molecular_surface_infile=blaster_files.thin_spheres_elec_molecular_surface_file,
                thin_spheres_outfile=blaster_files.thin_spheres_elec_file,
                distance_to_surface_parameter=flat_param_dict[
                    "thin_spheres_elec.distance_to_surface"
                ],
                penetration_parameter=flat_param_dict["thin_spheres_elec.penetration"],
            )
        )

        #
        steps.append(
            CloseSpheresGenerationStep(
                working_dir=working_dir,
                ligand_infile=blaster_files.ligand_hetatm_renamed_file,
                thin_spheres_infile=blaster_files.thin_spheres_elec_file,
                close_spheres_outfile=blaster_files.close_spheres_elec_file,
                distance_to_surface_parameter=flat_param_dict[
                    "thin_spheres_elec.distance_to_surface"
                ],
                penetration_parameter=flat_param_dict["thin_spheres_elec.penetration"],
                distance_to_ligand_parameter=flat_param_dict[
                    "thin_spheres_elec.distance_to_ligand"
                ],
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                working_dir=working_dir,
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
                working_dir=working_dir,
                charged_receptor_infile=blaster_files.charged_receptor_file,
                ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
                all_spheres_infile=blaster_files.all_spheres_file,
                low_dielectric_spheres_outfile=blaster_files.low_dielectric_spheres_file,
                **unpack_step_params(flat_param_dict, "low_dielectric_sphere_selection")
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                working_dir=working_dir,
                sph_infile=blaster_files.low_dielectric_spheres_file,
                pdb_outfile=blaster_files.low_dielectric_spheres_pdb_file,
            )
        )

        #
        spheres_pdb_file_for_electrostatics = (
            blaster_files.low_dielectric_spheres_pdb_file
        )

    #
    if flat_param_dict["thin_spheres_desolv.use"]:
        #
        steps.append(
            MolecularSurfaceGenerationStep(
                working_dir=working_dir,
                charged_receptor_infile=blaster_files.charged_receptor_file,
                binding_site_residues_infile=blaster_files.binding_site_residues_file,
                radii_infile=blaster_files.molecular_surface_radii_file,
                molecular_surface_outfile=blaster_files.thin_spheres_desolv_molecular_surface_file,
                molecular_surface_density_parameter=flat_param_dict["thin_spheres_desolv.molecular_surface_density"],
            )
        )

        #
        steps.append(
            ThinSpheresGenerationStep(
                working_dir=working_dir,
                molecular_surface_infile=blaster_files.thin_spheres_desolv_molecular_surface_file,
                thin_spheres_outfile=blaster_files.thin_spheres_desolv_file,
                distance_to_surface_parameter=flat_param_dict[
                    "thin_spheres_desolv.distance_to_surface"
                ],
                penetration_parameter=flat_param_dict[
                    "thin_spheres_desolv.penetration"
                ],
            )
        )

        #
        steps.append(
            CloseSpheresGenerationStep(
                working_dir=working_dir,
                ligand_infile=blaster_files.ligand_hetatm_renamed_file,
                thin_spheres_infile=blaster_files.thin_spheres_desolv_file,
                close_spheres_outfile=blaster_files.close_spheres_desolv_file,
                distance_to_surface_parameter=flat_param_dict[
                    "thin_spheres_desolv.distance_to_surface"
                ],
                penetration_parameter=flat_param_dict[
                    "thin_spheres_desolv.penetration"
                ],
                distance_to_ligand_parameter=flat_param_dict[
                    "thin_spheres_desolv.distance_to_ligand"
                ],
            )
        )

        #
        steps.append(
            SpheresToPDBConversionStep(
                working_dir=working_dir,
                sph_infile=blaster_files.close_spheres_desolv_file,
                pdb_outfile=blaster_files.close_spheres_desolv_pdb_file,
            )
        )

    #
    steps.append(
        BoxGenerationStep(
            working_dir=working_dir,
            charged_receptor_infile=blaster_files.charged_receptor_file,
            ligand_matching_spheres_infile=blaster_files.ligand_matching_spheres_file,
            box_outfile=blaster_files.box_file,
            **unpack_step_params(flat_param_dict, "box_generation")
        )
    )

    #
    steps.append(
        ReceptorTransformationForElectrostatics(
            working_dir=working_dir,
            charged_receptor_infile=blaster_files.charged_receptor_file,
            spheres_pdb_infile=spheres_pdb_file_for_electrostatics,
            receptor_low_dielectric_pdb_outfile=blaster_files.receptor_low_dielectric_pdb_file,
        )
    )

    #
    if flat_param_dict["thin_spheres_elec.use"]:
        steps.append(
            ElectrostaticsGridGenerationStepYesThinSpheres(
                working_dir=working_dir,
                receptor_low_dielectric_pdb_infile=blaster_files.receptor_low_dielectric_pdb_file,
                charge_infile=blaster_files.electrostatics_charge_file,
                radius_infile=blaster_files.electrostatics_radius_file,
                delphi_infile=blaster_files.electrostatics_delphi_file,
                box_infile=blaster_files.box_file,
                electrostatics_phi_outfile=blaster_files.electrostatics_phi_file,
                electrostatics_pdb_outfile=blaster_files.electrostatics_pdb_file,
                electrostatics_trim_phi_outfile=blaster_files.electrostatics_trim_phi_file,
                electrostatics_phi_size_outfile=blaster_files.electrostatics_phi_size_file,
                thin_spheres_elec_distance_to_surface_parameter=flat_param_dict[
                    "thin_spheres_elec.distance_to_surface"
                ],
                thin_spheres_elec_penetration_parameter=flat_param_dict[
                    "thin_spheres_elec.penetration"
                ],
                **unpack_step_params(flat_param_dict, "electrostatics_grid_gen")
            )
        )
    else:
        steps.append(
            ElectrostaticsGridGenerationStepNoThinSpheres(
                working_dir=working_dir,
                receptor_low_dielectric_pdb_infile=blaster_files.receptor_low_dielectric_pdb_file,
                charge_infile=blaster_files.electrostatics_charge_file,
                radius_infile=blaster_files.electrostatics_radius_file,
                delphi_infile=blaster_files.electrostatics_delphi_file,
                box_infile=blaster_files.box_file,
                electrostatics_phi_outfile=blaster_files.electrostatics_phi_file,
                electrostatics_pdb_outfile=blaster_files.electrostatics_pdb_file,
                electrostatics_trim_phi_outfile=blaster_files.electrostatics_trim_phi_file,
                electrostatics_phi_size_outfile=blaster_files.electrostatics_phi_size_file,
                **unpack_step_params(flat_param_dict, "electrostatics_grid_gen")
            )
        )

    #
    steps.append(
        VDWScoringGridGenerationStep(
            working_dir=working_dir,
            charged_receptor_infile=blaster_files.charged_receptor_file,
            vdw_parameters_infile=blaster_files.vdw_parameters_file,
            protein_table_infile=blaster_files.vdw_protein_table_file,
            box_infile=blaster_files.box_file,
            vdw_outfile=blaster_files.vdw_file,
            bump_map_outfile=blaster_files.vdw_bump_map_file,
            **unpack_step_params(flat_param_dict, "vdw_grid_gen")
        )
    )

    #
    if flat_param_dict["thin_spheres_desolv.use"]:
        steps.append(
            ReceptorTransformationForLigandDesolvationYesThinSpheres(
                working_dir=working_dir,
                charged_receptor_pdb_infile=blaster_files.charged_receptor_file,
                close_spheres_desolv_pdb_infile=blaster_files.close_spheres_desolv_pdb_file,
                charged_receptor_desolv_pdb_outfile=blaster_files.charged_receptor_desolv_pdb_file,
            )
        )
    else:
        steps.append(
            ReceptorTransformationForLigandDesolvationNoThinSpheres(
                working_dir=working_dir,
                charged_receptor_pdb_infile=blaster_files.charged_receptor_file,
                charged_receptor_desolv_pdb_outfile=blaster_files.charged_receptor_desolv_pdb_file,
            )
        )

    #
    steps.append(
        HydrogenAtomLigandDesolvationScoringGridGenerationStep(
            working_dir=working_dir,
            box_infile=blaster_files.box_file,
            receptor_pdb_infile=blaster_files.charged_receptor_desolv_pdb_file,
            ligand_desolvation_outfile=blaster_files.ligand_desolvation_hydrogen_file,
            thin_spheres_desolv_use_parameter=flat_param_dict[
                "thin_spheres_desolv.use"
            ],
            thin_spheres_desolv_distance_to_surface_parameter=flat_param_dict[
                "thin_spheres_desolv.distance_to_surface"
            ],
            thin_spheres_desolv_penetration_parameter=flat_param_dict[
                "thin_spheres_desolv.penetration"
            ],
            **unpack_step_params(flat_param_dict, "desolv_grid_gen")
        )
    )

    #
    steps.append(
        HeavyAtomLigandDesolvationScoringGridGenerationStep(
            working_dir=working_dir,
            box_infile=blaster_files.box_file,
            receptor_pdb_infile=blaster_files.charged_receptor_desolv_pdb_file,
            ligand_desolvation_outfile=blaster_files.ligand_desolvation_heavy_file,
            thin_spheres_desolv_use_parameter=flat_param_dict[
                "thin_spheres_desolv.use"
            ],
            thin_spheres_desolv_distance_to_surface_parameter=flat_param_dict[
                "thin_spheres_desolv.distance_to_surface"
            ],
            thin_spheres_desolv_penetration_parameter=flat_param_dict[
                "thin_spheres_desolv.penetration"
            ],
            **unpack_step_params(flat_param_dict, "desolv_grid_gen")
        )
    )

    steps.append(
        VisualizationStep(
            working_dir=working_dir,
            vdw_infile=blaster_files.vdw_file,
            vdw_bump_map_infile=blaster_files.vdw_bump_map_file,
            lig_desolv_heavy_infile=blaster_files.ligand_desolvation_heavy_file,
            trim_electrostatics_phi_infile=blaster_files.electrostatics_trim_phi_file,
            matching_spheres_infile=blaster_files.matching_spheres_file,
            vdw_repulsive_dx_outfile=blaster_files.vdw_repulsive_dx_file,
            vdw_attractive_dx_outfile=blaster_files.vdw_attractive_dx_file,
            vdw_dx_outfile=blaster_files.vdw_dx_file,
            lig_desolv_dx_outfile=blaster_files.ligand_desolvation_dx_file,
            trim_electrostatics_dx_outfile=blaster_files.trim_electrostatics_dx_file,
            matching_spheres_outfile=blaster_files.matching_spheres_pdb_file
        )
    )

    return tuple(steps)


def load_steps(job_dir_path, config_file_path):
    """Returns the blaster files, flattened config parameters, and blaster steps of a job."""
    working_dir = WorkingDir(
        path=os.path.join(job_dir_path, Blastermaster.WORKING_DIR_NAME),
        create=False,
        reset=False,
    )
    blaster_files = BlasterFiles(working_dir=working_dir)
    config = BlastermasterParametersConfiguration(config_file_path)
    flat_param_dict = flatten_and_parameter_cast_param_dict(config.param_dict)
    dock_files_generation_flat_param_dict = {key.replace('dock_files_generation.', ''): value for key, value in flat_param_dict.items() if key.startswith('dock_files_generation.')}
    steps = get_blaster_steps(blaster_files, dock_files_generation_flat_param_dict, working_dir)
    return blaster_files, flat_param_dict, steps


def convert_line_endings_to_unix(file_path):
    """Programs compiled for Windows write CRLF line endings; DOCK expects LF."""
    with open(file_path, "rb") as f:
        content = f.read()
    if b"\r\n" in content:
        with open(file_path, "wb") as f:
            f.write(content.replace(b"\r\n", b"\n"))


class Blastermaster(Script):
    JOB_DIR_NAME = "blastermaster_job"
    CONFIG_FILE_NAME = "blastermaster_config.yaml"
    DEFAULT_CONFIG_FILE_PATH = os.path.join(
        os.path.dirname(BLASTERMASTER_INIT_FILE_PATH),
        "default_blastermaster_config.yaml",
    )
    WORKING_DIR_NAME = "working"
    DOCK_FILES_DIR_NAME = "dockfiles"
    VISUALIZATION_FILES_DIR_NAME = "visualization"
    DEFAULT_FILES_DIR_PATH = os.path.dirname(DEFAULTS_INIT_FILE_PATH)

    def __init__(self):
        super().__init__()

    def new(self, job_dir_path=JOB_DIR_NAME, overwrite=False):
        # create job dir
        job_dir = Dir(path=job_dir_path, create=True, reset=False)

        # create working dir & copy in blaster files
        blaster_file_names = list(BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT.values())
        backup_blaster_file_paths = [
            os.path.join(self.DEFAULT_FILES_DIR_PATH, blaster_file_name)
            for blaster_file_name in blaster_file_names
        ]
        blaster_file_names_in_cwd = [f for f in blaster_file_names if os.path.isfile(f)]
        files_to_copy_str = "\n\t".join(blaster_file_names_in_cwd)
        if blaster_file_names_in_cwd:
            logger.info(
                f"Copying the following files from current directory into job working directory:\n\t{files_to_copy_str}"
            )
        working_dir = WorkingDir(
            path=os.path.join(job_dir.path, self.WORKING_DIR_NAME),
            create=True,
            reset=False,
            files_to_copy_in=blaster_file_names_in_cwd,
            new_file_names=blaster_file_names_in_cwd,
            backup_files_to_copy_in=backup_blaster_file_paths,
            new_backup_file_names=blaster_file_names,
        )

        # create dock files dir
        dock_files_dir = Dir(
            path=os.path.join(job_dir.path, self.DOCK_FILES_DIR_NAME),
            create=True,
            reset=False,
        )
        visualization_dir = Dir(
            path=os.path.join(job_dir_path, self.VISUALIZATION_FILES_DIR_NAME),
            create=True,
            reset=False,
        )

        # write fresh config file from default file
        save_path = os.path.join(job_dir.path, self.CONFIG_FILE_NAME)
        BlastermasterParametersConfiguration.write_config_file(
            save_path, self.DEFAULT_CONFIG_FILE_PATH, overwrite=overwrite
        )

    def run(
        self,
        job_dir_path=".",
        config_file_path=None,
        #use_graph_state=True,  # TODO
        #write_graph_image=False,  # TODO
    ):
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
        dock_files_dir = Dir(
            path=os.path.join(job_dir.path, self.DOCK_FILES_DIR_NAME),
            create=True,
            reset=True,
        )  # reset dock files dir in case re-running
        visualization_dir = Dir(
            path=os.path.join(job_dir_path, self.VISUALIZATION_FILES_DIR_NAME),
            create=True,
            reset=True,
        )

        #
        indock_file = IndockFile(
            path=os.path.join(dock_files_dir.path, INDOCK_FILE_NAME)
        )

        # load config file & blaster steps
        logger.info("Loading config file")
        blaster_files, flat_param_dict, steps = load_steps(job_dir.path, config_file_path)

        # get params as str
        config_params_str = "\n".join(
            [
                f"{param_name}: {param.value}"
                for param_name, param in flat_param_dict.items()
            ]
        )
        logger.info(f"Parameters:\n{config_params_str}")

        # reset step dirs
        for step in steps:
            if not step.is_done:
                step.step_dir.delete()

        # run steps
        logger.info("Running blaster steps in sequence.")
        for step in steps:
            step.run()

        # copy dock files to dock files directory
        logger.info("Copying dock files to dock files directory")
        for dock_file in astuple(blaster_files.dock_files):
            dock_file_path = os.path.join(dock_files_dir.path, dock_file.name)
            File.copy_file(dock_file.path, dock_file_path)
            if dock_file.identifier not in BINARY_DOCK_FILE_IDENTIFIERS:
                convert_line_endings_to_unix(dock_file_path)

        # copy visualization files to visualization directory
        logger.info("Copying visualization files to visualization directory")
        for visualization_file in astuple(blaster_files.visualization_files):
            File.copy_file(
                visualization_file.path, os.path.join(visualization_dir.path, visualization_file.name)
            )

        # write INDOCK file
        logger.info("Making indock file")
        indock_file_generation_flat_param_dict = {key.replace('indock_file_generation.', ''): value for key, value in flat_param_dict.items() if key.startswith('indock_file_generation.')}
        indock_file.write(
            blaster_files.dock_files, indock_file_generation_flat_param_dict, dock_files_dir.name
        )

