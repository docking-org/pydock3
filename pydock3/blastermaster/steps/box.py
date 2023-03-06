import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BoxGenerationStep(BlasterStep):

    MARGIN = 10.0

    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        box_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
                (ligand_matching_spheres_infile, "ligand_matching_spheres_infile", None),
            ],
            outfile_tuples=[
                (box_outfile, "box_outfile", None),
            ],
            parameter_tuples=[],
            program_file_path=ProgramFilePaths.MAKEBOX_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run the makebox.smallokay.pl perl script to make box surrounding binding
        site"""

        run_str = f"{self.program_file.path} {self.infiles.ligand_matching_spheres_infile.name} {self.infiles.charged_receptor_infile.name} {self.outfiles.box_outfile.name} {self.MARGIN}"
        self.run_command(run_str)
