#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab
# Trent E. Balius modified Sept 2013.

import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorMostOccupiedResiduesRenamingStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        receptor_infile,
        receptor_most_occupied_residues_renamed_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (receptor_infile, "receptor_infile"),
        )

        #
        self.process_outfiles(
            (
                receptor_most_occupied_residues_renamed_outfile,
                "receptor_most_occupied_residues_renamed_outfile",
            ),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        #
        pdb.move_columns(
            input_pdb_file_path=self.infiles.receptor_infile.path,
            output_pdb_file_path=self.outfiles.receptor_most_occupied_residues_renamed_outfile.path,
        )
        logger.info("Getting the most occupied position of each residue and renaming")
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
