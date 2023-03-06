import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BindingSiteResiduesSelectionStep(BlasterStep):
    class MandatoryFileNames:
        RECEPTOR_FILE_NAME = "rec.pdb"
        LIGAND_FILE_NAME = "xtal-lig.pdb"
        BINDING_SITE_RESIDUES_FILE_NAME = "rec.site"

    def __init__(
        self,
        working_dir,
        receptor_infile,
        ligand_infile,
        filt_parameters_infile,
        binding_site_residues_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_infile, "receptor_infile", self.MandatoryFileNames.RECEPTOR_FILE_NAME),
                (ligand_infile, "ligand_infile", self.MandatoryFileNames.LIGAND_FILE_NAME),
                (filt_parameters_infile, "filt_parameters_infile", None),
            ],
            outfile_tuples=[
                (binding_site_residues_outfile, "binding_site_residues_outfile", self.MandatoryFileNames.BINDING_SITE_RESIDUES_FILE_NAME),
            ],
            parameter_tuples=[],
            program_file_path=ProgramFilePaths.FILT_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """just run the filt.exe program to produce list of binding site residues
        The program filt.exe presumably outputs a file called rec.site in the current directory.
        """
        #
        run_str = (
            f"{self.program_file.path} < {self.infiles.filt_parameters_infile.name}"
        )
        self.run_command(run_str)
