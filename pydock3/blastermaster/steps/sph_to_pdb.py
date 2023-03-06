import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class SpheresToPDBConversionStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        sph_infile,
        pdb_outfile,
    ):
        #
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (sph_infile, "sph_infile", None),
            ],
            outfile_tuples=[
                (pdb_outfile, "pdb_outfile", None),
            ],
            parameter_tuples=[],
            program_file_path=ProgramFilePaths.DOSHOWSPH_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """convert spheres to pdb file"""

        run_str = f"{self.program_file.path} {self.infiles.sph_infile.name} 1 {self.outfiles.pdb_outfile.name}"
        self.run_command(run_str)
