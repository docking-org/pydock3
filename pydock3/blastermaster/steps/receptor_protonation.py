import logging

import yaml

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import File
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorProtonationStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        receptor_infile,
        add_h_dict_infile,
        residue_code_polar_h_yaml_infile,
        charged_receptor_outfile,
        reduce_options_parameter,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_infile, "receptor_infile", None),
                (add_h_dict_infile, "add_h_dict_infile", None),
                (residue_code_polar_h_yaml_infile, "residue_code_polar_h_yaml_infile", None),
            ],
            outfile_tuples=[
                (charged_receptor_outfile, "charged_receptor_outfile", None),
            ],
            parameter_tuples=[
                (reduce_options_parameter, "reduce_options_parameter"),
            ],
            program_file_path=ProgramFilePaths.REDUCE_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run REDUCE to produce a pdb with hydrogens.
        Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.
        then run script to remove nonpolar hydrogens & rename"""
        #
        charged_receptor_full_h_file_path = (
            f"{self.outfiles.charged_receptor_outfile.path}.fullh"
        )
        charged_receptor_full_h_file_name = File.get_file_name_of_file(
            charged_receptor_full_h_file_path
        )
        run_str = f"{self.program_file.path} -db {self.infiles.add_h_dict_infile.name} {self.parameters.reduce_options_parameter.value} {self.infiles.receptor_infile.name} > {charged_receptor_full_h_file_name}"
        self.run_command(run_str)

        # remove extraneous output from charged_receptor_full_h_file_path
        run_str = f"sed -i 's/\s*new\s*//g' {charged_receptor_full_h_file_name} ; sed -i '/^USER.*/d' {charged_receptor_full_h_file_name}"
        self.run_command(run_str)

        # remove nonpolar hydrogens
        pdb_d = pdb.PDBData(charged_receptor_full_h_file_path, ignore_waters=False)
        with open(self.infiles.residue_code_polar_h_yaml_infile.path, 'r') as f:
            residue_code_to_polar_hydrogens_dict = yaml.safe_load(f)
            pdb_d.remove_apolar_hydrogen(residue_code_to_polar_hydrogens_dict)

        #
        charged_receptor_polar_h_file_path = (
            f"{self.outfiles.charged_receptor_outfile.path}.polarH"
        )
        pdb_d.write(charged_receptor_polar_h_file_path)

        # rename histidines and cysteines
        pdb_d = pdb.PDBData(charged_receptor_polar_h_file_path, ignore_waters=False)
        pdb_d.rename_histidines()
        pdb_d.rename_cysteines()

        #
        pdb_d.write(self.outfiles.charged_receptor_outfile.path)
