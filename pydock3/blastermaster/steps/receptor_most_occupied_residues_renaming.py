import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorMostOccupiedResiduesRenamingStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        receptor_infile,
        receptor_most_occupied_residues_renamed_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_infile, "receptor_infile", None),
            ],
            outfile_tuples=[
                (receptor_most_occupied_residues_renamed_outfile, "receptor_most_occupied_residues_renamed_outfile", None),
            ],
            parameter_tuples=[],
            program_file_path=None,
        )

    @BlasterStep.handle_run_func
    def run(self):
        #
        pdb.move_columns(
            input_pdb_file_path=self.infiles.receptor_infile.path,
            output_pdb_file_path=self.outfiles.receptor_most_occupied_residues_renamed_outfile.path,
        )
        logger.debug("Getting the most occupied position of each residue and renaming")
        pdb.most_occupied(
            pdb_file_path=self.outfiles.receptor_most_occupied_residues_renamed_outfile.path,
            outfile_path=self.outfiles.receptor_most_occupied_residues_renamed_outfile.path,
        )
        pdb_pre_rename_residues = pdb.PDBData(
            self.outfiles.receptor_most_occupied_residues_renamed_outfile.path
        )
        pdb_pre_rename_residues.replace_alt_chars(" ")
        pdb_pre_rename_residues.delete_insertion_codes()
        pdb_pre_rename_residues.fix_chain_ids()
        pdb_pre_rename_residues.write(
            self.outfiles.receptor_most_occupied_residues_renamed_outfile.path
        )
