# Ryan G. Coleman

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile

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
        step_dir,
        receptor_infile,
        ligand_infile,
        filt_parameters_infile,
        binding_site_residues_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(path=ProgramFilePaths.FILT_PROGRAM_FILE_PATH)

        #
        self.process_infiles(
            (receptor_infile, "receptor_infile"),
            (ligand_infile, "ligand_infile"),
            (filt_parameters_infile, "filt_parameters_infile"),
            new_file_names_tuple=(self.MandatoryFileNames.RECEPTOR_FILE_NAME, self.MandatoryFileNames.LIGAND_FILE_NAME, filt_parameters_infile.name),
        )

        #
        self.process_outfiles(
            (binding_site_residues_outfile, "binding_site_residues_outfile"),
            new_file_names_tuple=(self.MandatoryFileNames.BINDING_SITE_RESIDUES_FILE_NAME,),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        """just run the filt.exe program to produce list of binding site residues
        The program filt.exe presumably outputs a file called rec.site in the current directory.
        """
        #
        run_str = f"{self.program_file.path} < {self.infiles.filt_parameters_infile.name}"
        self.run_command(run_str)
