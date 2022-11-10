# Ryan G. Coleman, Brian K. Shoichet Lab
# Trent E. Balius modified Sept 2013.

import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile, File
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorProtonationStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        receptor_infile,
        add_h_dict_infile,
        charged_receptor_outfile,
        reduce_options_parameter,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(path=ProgramFilePaths.REDUCE_PROGRAM_FILE_PATH)

        #
        self.process_infiles(
            (receptor_infile, "receptor_infile"),
            (add_h_dict_infile, "add_h_dict_infile"),
        )

        #
        self.process_outfiles(
            (charged_receptor_outfile, "charged_receptor_outfile"),
        )

        #
        self.process_parameters(
            (reduce_options_parameter, "reduce_options_parameter"),
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run REDUCE to produce a pdb with hydrogens.
        Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.
        then run script to remove nonpolar hydrogens & rename"""
        #
        charged_receptor_full_h_file_path = f"{self.outfiles.charged_receptor_outfile.path}.fullh"
        charged_receptor_full_h_file_name = File.get_file_name_of_file(charged_receptor_full_h_file_path)
        run_str = f"{self.program_file.path} -db {self.infiles.add_h_dict_infile.name} {self.parameters.reduce_options_parameter.value} {self.infiles.receptor_infile.name} > {charged_receptor_full_h_file_name}"
        self.run_command(run_str)

        # remove extraneous output from charged_receptor_full_h_file_path
        run_str = f"sed -i 's/\s*new\s*//g' {charged_receptor_full_h_file_name} ; sed -i '/^USER.*/d' {charged_receptor_full_h_file_name}"
        self.run_command(run_str)

        # remove nonpolar hydrogens
        pdb_d = pdb.PDBData(charged_receptor_full_h_file_path, ignore_waters=False)
        pdb_d.remove_apolar_hydrogen()

        #
        charged_receptor_polar_h_file_path = f"{self.outfiles.charged_receptor_outfile.path}.polarH"
        pdb_d.write(charged_receptor_polar_h_file_path)

        # rename histidines and cysteines
        pdb_d = pdb.PDBData(charged_receptor_polar_h_file_path, ignore_waters=False)
        pdb_d.rename_histidines()
        pdb_d.rename_cysteines()

        #
        pdb_d.write(self.outfiles.charged_receptor_outfile.path)
