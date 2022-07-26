#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from ucsfdock.blastermaster.util import ProgramFilePaths, BlasterStep
from ucsfdock.files import ProgramFile


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandPDBToSpheresConversionStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        pdb_infile,
        sph_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.PDBTOSPH_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (pdb_infile, "pdb_infile"),
        )

        #
        self.process_outfiles(
            (sph_outfile, "sph_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        """run the pdbtosph program to turn the ligand into spheres"""
        #
        run_str = f"{self.program_file.path} {self.infiles.pdb_infile.name} {self.outfiles.sph_outfile.name}"
        self.run_command(run_str)
