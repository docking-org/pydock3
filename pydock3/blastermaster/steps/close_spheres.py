#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class CloseSpheresGenerationStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        ligand_infile,
        thin_spheres_infile,
        close_spheres_outfile,
        distance_to_surface_parameter,
        penetration_parameter,
        distance_to_ligand_parameter,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.CLOSE_SPH_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (ligand_infile, "ligand_infile"),
            (thin_spheres_infile, "thin_spheres_infile"),
        )

        #
        self.process_outfiles(
            (close_spheres_outfile, "close_spheres_outfile"),
        )

        #
        self.process_parameters(
            (distance_to_surface_parameter, "distance_to_surface_parameter,"),
            (penetration_parameter, "penetration_parameter,"),
            (distance_to_ligand_parameter, "distance_to_ligand_parameter,"),
        )

    @BlasterStep.handle_run_func
    def run(self):
        #
        run_str = f"python {self.program_file.path} {self.infiles.thin_spheres_infile.name} {self.infiles.ligand_infile.name} {self.outfiles.close_spheres_outfile.name} {self.parameters.distance_to_ligand_parameter.value} {self.parameters.distance_to_surface_parameter.value + self.parameters.penetration_parameter.value}"
        self.run_command(run_str)
