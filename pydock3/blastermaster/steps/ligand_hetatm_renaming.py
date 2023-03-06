import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandHetatmRenamingStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        ligand_infile,
        ligand_hetatm_renamed_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (ligand_infile, "ligand_infile", None),
            ],
            outfile_tuples=[
                (ligand_hetatm_renamed_outfile, "ligand_hetatm_renamed_outfile", None),
            ],
            parameter_tuples=[],
            program_file_path=None,
        )

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
