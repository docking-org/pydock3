#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from ucsfdock.blastermaster.util import ProgramFilePaths, BlasterStep
from ucsfdock.files import ProgramFile


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class SpheresToPDBConversionStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        sph_infile,
        pdb_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.DOSHOWSPH_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (sph_infile, "sph_infile"),
        )

        #
        self.process_outfiles(
            (pdb_outfile, "pdb_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        """convert spheres to pdb file"""

        run_str = f"{self.program_file.path} {self.infiles.sph_infile.name} 1 {self.outfiles.pdb_outfile.name}"
        self.run_command(run_str)
