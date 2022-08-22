#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ThinSpheresGenerationStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        molecular_surface_infile,
        thin_spheres_outfile,
        distance_to_surface_parameter,
        penetration_parameter,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.THIN_SPHERES_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (molecular_surface_infile, "molecular_surface_infile"),
        )

        #
        self.process_outfiles(
            (thin_spheres_outfile, "thin_spheres_outfile"),
        )

        #
        self.process_parameters(
            (distance_to_surface_parameter, "distance_to_surface_parameter"),
            (penetration_parameter, "penetration_parameter"),
        )

    @BlasterStep.handle_run_func
    def run(self):
        run_str = f"{self.program_file.path} -i {self.infiles.molecular_surface_infile.name} -o {self.outfiles.thin_spheres_outfile.name} -d {self.parameters.distance_to_surface_parameter.value} -s {self.parameters.distance_to_surface_parameter.value + self.parameters.penetration_parameter.value}"
        self.run_command(run_str)
