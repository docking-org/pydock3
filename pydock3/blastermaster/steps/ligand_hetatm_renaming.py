# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandHetatmRenamingStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        ligand_infile,
        ligand_hetatm_renamed_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (ligand_infile, "ligand_infile"),
        )

        #
        self.process_outfiles(
            (ligand_hetatm_renamed_outfile, "ligand_hetatm_renamed_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        #
        logger.info("Renaming HETATM to ATOM")
        pdb.move_columns(
            input_pdb_file_path=self.infiles.ligand_infile.path,
            output_pdb_file_path=self.outfiles.ligand_hetatm_renamed_outfile.path,
        )
        pdb_ligand = pdb.PDBData(self.outfiles.ligand_hetatm_renamed_outfile.path)
        pdb_ligand.replace_hetatm_with_atom()
        pdb_ligand.write(self.outfiles.ligand_hetatm_renamed_outfile.path)
