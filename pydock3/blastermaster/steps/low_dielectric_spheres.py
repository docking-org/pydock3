# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LowDielectricSpheresSelectionStep(BlasterStep):

    MIN_NUM_SPHERES = 25

    def __init__(
        self,
        step_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        all_spheres_infile,
        low_dielectric_spheres_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.MAKESPHERES1_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (charged_receptor_infile, "charged_receptor_infile"),
            (ligand_matching_spheres_infile, "ligand_matching_spheres_infile"),
            (all_spheres_infile, "all_spheres_infile"),
        )

        #
        self.process_outfiles(
            (low_dielectric_spheres_outfile, "low_dielectric_spheres_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        """run the makespheres1.cli.pl perl script to make low dielectric spheres"""

        run_str = f"{self.program_file.path} {self.infiles.ligand_matching_spheres_infile.name} {self.infiles.all_spheres_infile.name} {self.infiles.charged_receptor_infile.name} {self.outfiles.low_dielectric_spheres_outfile.name} {self.MIN_NUM_SPHERES}"
        self.run_command(run_str)
