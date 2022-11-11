# Ryan G. Coleman, Brian K. Shoichet Lab

import os
import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile, File

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

    ATOM_TYPE_TO_RADIUS_DICT = {
        AtomTypes.HYDROGEN: 1.0,
        AtomTypes.HEAVY: 1.8,
    }
    PROBE_RADIUS = 1.4  # radius of water

    def __init__(
        self,
        step_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        other_radius_parameter,
        atom_type,
    ):
        #
        if atom_type not in self.ATOM_TYPE_TO_RADIUS_DICT:
            logger.exception(f"atom_type must be one of: {list(self.ATOM_TYPE_TO_RADIUS_DICT.keys())}")
            raise

        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.SOLVMAP_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (receptor_pdb_infile, "receptor_pdb_infile"),
            (box_infile, "box_infile"),
        )

        #
        self.process_outfiles(
            (ligand_desolvation_outfile, "ligand_desolvation_outfile"),
        )

        #
        self.process_parameters(
            (thin_spheres_desolv_use_parameter, "thin_spheres_desolv_use_parameter"),
            (thin_spheres_desolv_distance_to_surface_parameter, "thin_spheres_desolv_distance_to_surface_parameter"),
            (thin_spheres_desolv_penetration_parameter, "thin_spheres_desolv_penetration_parameter"),
            (other_radius_parameter, "other_radius_parameter"),
        )

        # misc.
        self.atom_type = atom_type

    @BlasterStep.handle_run_func
    def run(self):
        """run the solvmap program"""

        # make solvmap parameters file
        solvmap_parameters_file = File(path=os.path.join(self.step_dir.path, self.MandatoryFileNames.SOLVMAP_PARAMETERS_FILE_NAME))
        with open(solvmap_parameters_file.path, "w") as f:
            f.write(f"{self.infiles.receptor_pdb_infile.name}\n")  # receptor file name
            f.write(f"{self.outfiles.ligand_desolvation_outfile.name}\n")  # output file name
            if self.parameters.thin_spheres_desolv_use_parameter.value:
                other_radius = self.parameters.thin_spheres_desolv_distance_to_surface_parameter.value + self.parameters.thin_spheres_desolv_penetration_parameter.value
            else:
                other_radius = self.parameters.other_radius_parameter.value
            f.write(
                "1.60,1.65,1.90,1.90,1.90,%3.2f\n" % other_radius
            )  # radius of O,N,C,S,P,other
            f.write(f"{self.PROBE_RADIUS}\n")  # probe radius
            f.write("2\n")  # grid resolution
            f.write(f"{self.infiles.box_infile.name}\n")  # box file, extent of grids
            f.write(f"{self.ATOM_TYPE_TO_RADIUS_DICT[self.atom_type]}\n")  # born radius

        #
        self.log_parameters_file(solvmap_parameters_file)

        # run
        run_str = f"{self.program_file.path}"
        self.run_command(run_str)


class HydrogenAtomLigandDesolvationScoringGridGenerationStep(LigandDesolvationScoringGridGenerationStep):
    def __init__(
        self,
        step_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        other_radius_parameter,
    ):
        super().__init__(
            step_dir=step_dir,
            receptor_pdb_infile=receptor_pdb_infile,
            box_infile=box_infile,
            ligand_desolvation_outfile=ligand_desolvation_outfile,
            thin_spheres_desolv_use_parameter=thin_spheres_desolv_use_parameter,
            thin_spheres_desolv_distance_to_surface_parameter=thin_spheres_desolv_distance_to_surface_parameter,
            thin_spheres_desolv_penetration_parameter=thin_spheres_desolv_penetration_parameter,
            other_radius_parameter=other_radius_parameter,
            atom_type=super().AtomTypes.HYDROGEN
        )


class HeavyAtomLigandDesolvationScoringGridGenerationStep(LigandDesolvationScoringGridGenerationStep):
    def __init__(
        self,
        step_dir,
        receptor_pdb_infile,
        box_infile,
        ligand_desolvation_outfile,
        thin_spheres_desolv_use_parameter,
        thin_spheres_desolv_distance_to_surface_parameter,
        thin_spheres_desolv_penetration_parameter,
        other_radius_parameter,
    ):
        super().__init__(
            step_dir=step_dir,
            receptor_pdb_infile=receptor_pdb_infile,
            box_infile=box_infile,
            ligand_desolvation_outfile=ligand_desolvation_outfile,
            thin_spheres_desolv_use_parameter=thin_spheres_desolv_use_parameter,
            thin_spheres_desolv_distance_to_surface_parameter=thin_spheres_desolv_distance_to_surface_parameter,
            thin_spheres_desolv_penetration_parameter=thin_spheres_desolv_penetration_parameter,
            other_radius_parameter=other_radius_parameter,
            atom_type=super().AtomTypes.HEAVY
        )
