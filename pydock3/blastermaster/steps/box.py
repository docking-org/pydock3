# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BoxGenerationStep(BlasterStep):

    MARGIN = 10.0

    def __init__(
        self,
        step_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        box_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(path=ProgramFilePaths.MAKEBOX_PROGRAM_FILE_PATH)

        #
        self.process_infiles(
            (charged_receptor_infile, "charged_receptor_infile"),
            (ligand_matching_spheres_infile, "ligand_matching_spheres_infile"),
        )

        #
        self.process_outfiles(
            (box_outfile, "box_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        """run the makebox.smallokay.pl perl script to make box surrounding binding
        site"""

        run_str = f"{self.program_file.path} {self.infiles.ligand_matching_spheres_infile.name} {self.infiles.charged_receptor_infile.name} {self.outfiles.box_outfile.name} {self.MARGIN}"
        self.run_command(run_str)
