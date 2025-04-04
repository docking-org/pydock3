import os
import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile, File
from pydock3.config import Parameter

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandDesolvationScoringGridGenerationStep(BlasterStep):
    class MandatoryFileNames:
        SOLVMAP_PARAMETERS_FILE_NAME = (
            "INSEV"  # hardwired, stupid. fix in solvmap sometime
        )

    class AtomTypes:
        HYDROGEN = 1
        HEAVY = 2

    # Defaults that can be overwritten in the config
    OTHER_RADIUS = Parameter("dock_files_generation.desolv_grid_gen.other_radius", 1.0)
    PROBE_RADIUS = Parameter("dock_files_generation.desolv_grid_gen.probe_radius", 1.4)  # radius of water
    HYDROGEN_RADIUS = Parameter("dock_files_generation.desolv_grid_gen.hydrogen_radius", 1.0) 
    HEAVY_RADIUS = Parameter("dock_files_generation.desolv_grid_gen.heavy_radius", 1.8)
    SUBMIT_TO_SCHEDULER = Parameter("dock_files_generation.desolv_grid_gen.submit_to_scheduler", False)

    def __init__(
        self,
        working_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        atom_type,
        other_radius_parameter=OTHER_RADIUS,
        probe_radius_parameter=PROBE_RADIUS,
        hydrogen_radius_parameter=HYDROGEN_RADIUS,
        heavy_radius_parameter=HEAVY_RADIUS,
        submit_to_scheduler_parameter=SUBMIT_TO_SCHEDULER,
    ):

        self.ATOM_TYPE_TO_RADIUS_DICT = {
            self.AtomTypes.HYDROGEN: hydrogen_radius_parameter.value,
            self.AtomTypes.HEAVY: heavy_radius_parameter.value 
        }

        #
        if atom_type not in self.ATOM_TYPE_TO_RADIUS_DICT:
            logger.exception(
                f"atom_type must be one of: {list(self.ATOM_TYPE_TO_RADIUS_DICT.keys())}"
            )
            raise

        #
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_pdb_infile, "receptor_pdb_infile", None),
                (box_infile, "box_infile", None),
            ],
            outfile_tuples=[
                (ligand_desolvation_outfile, "ligand_desolvation_outfile", None),
            ],
            parameter_tuples=[
                (thin_spheres_desolv_use_parameter, "thin_spheres_desolv_use_parameter"),
                (thin_spheres_desolv_distance_to_surface_parameter, "thin_spheres_desolv_distance_to_surface_parameter"),
                (thin_spheres_desolv_penetration_parameter, "thin_spheres_desolv_penetration_parameter"),
                (other_radius_parameter, "other_radius_parameter"),
                (probe_radius_parameter, "probe_radius_parameter"),
                (hydrogen_radius_parameter, "hydrogen_radius_parameter"),
                (heavy_radius_parameter, "heavy_radius_parameter"),

            ],
            program_file_path=ProgramFilePaths.SOLVMAP_PROGRAM_FILE_PATH,
            dockopt_submit_to_scheduler=submit_to_scheduler_parameter.value,
        )

        # misc.
        self.atom_type = atom_type

    @BlasterStep.handle_run_func
    def run(self):
        """run the solvmap program"""

        # make solvmap parameters file
        solvmap_parameters_file = File(
            path=os.path.join(
                self.step_dir.path, self.MandatoryFileNames.SOLVMAP_PARAMETERS_FILE_NAME
            )
        )
        with open(solvmap_parameters_file.path, "w") as f:
            f.write(f"{self.infiles.receptor_pdb_infile.name}\n")  # receptor file name
            f.write(
                f"{self.outfiles.ligand_desolvation_outfile.name}\n"
            )  # output file name
            if self.parameters.thin_spheres_desolv_use_parameter.value:
                other_radius = (
                    self.parameters.thin_spheres_desolv_distance_to_surface_parameter.value
                    + self.parameters.thin_spheres_desolv_penetration_parameter.value
                )
            else:
                other_radius = self.parameters.other_radius_parameter.value
            f.write(
                "1.60,1.65,1.90,1.90,1.90,%3.2f\n" % other_radius
            )  # radius of O,N,C,S,P,other
            f.write(f"{self.parameters.probe_radius_parameter.value}\n")  # probe radius
            f.write("2\n")  # grid resolution
            f.write(f"{self.infiles.box_infile.name}\n")  # box file, extent of grids
            f.write(f"{self.ATOM_TYPE_TO_RADIUS_DICT[self.atom_type]}\n")  # born radius

        #
        self.log_parameters_file(solvmap_parameters_file)

        # run
        run_str = f"{self.program_file.path}"
        self.run_command(run_str)


class HydrogenAtomLigandDesolvationScoringGridGenerationStep(
    LigandDesolvationScoringGridGenerationStep
):
    def __init__(
        self,
        working_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        **kwargs
    ):
        super().__init__(
            working_dir=working_dir,
            receptor_pdb_infile=receptor_pdb_infile,
            box_infile=box_infile,
            ligand_desolvation_outfile=ligand_desolvation_outfile,
            thin_spheres_desolv_use_parameter=thin_spheres_desolv_use_parameter,
            thin_spheres_desolv_distance_to_surface_parameter=thin_spheres_desolv_distance_to_surface_parameter,
            thin_spheres_desolv_penetration_parameter=thin_spheres_desolv_penetration_parameter,
            atom_type=super().AtomTypes.HYDROGEN,
            **kwargs
        )


class HeavyAtomLigandDesolvationScoringGridGenerationStep(
    LigandDesolvationScoringGridGenerationStep
):
    def __init__(
        self,
        working_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        **kwargs
    ):
        super().__init__(
            working_dir=working_dir,
            receptor_pdb_infile=receptor_pdb_infile,
            box_infile=box_infile,
            ligand_desolvation_outfile=ligand_desolvation_outfile,
            thin_spheres_desolv_use_parameter=thin_spheres_desolv_use_parameter,
            thin_spheres_desolv_distance_to_surface_parameter=thin_spheres_desolv_distance_to_surface_parameter,
            thin_spheres_desolv_penetration_parameter=thin_spheres_desolv_penetration_parameter,
            atom_type=super().AtomTypes.HEAVY,
            **kwargs
        )
