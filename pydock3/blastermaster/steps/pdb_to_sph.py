import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandPDBToSpheresConversionStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        pdb_infile,
        sph_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (pdb_infile, "pdb_infile", None),
            ],
            outfile_tuples=[
                (sph_outfile, "sph_outfile", None),
            ],
            parameter_tuples=[],
            program_file_path=ProgramFilePaths.PDBTOSPH_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run the pdbtosph program to turn the ligand into spheres"""
        #
        run_str = f"{self.program_file.path} {self.infiles.pdb_infile.name} {self.outfiles.sph_outfile.name}"
        self.run_command(run_str)
