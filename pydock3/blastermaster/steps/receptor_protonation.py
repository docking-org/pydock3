import logging
import re
import shlex

import yaml

from pydock3.blastermaster.util import program_path, BlasterStep
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
        self.run_program(
            [
                program_path("reduce"),
                "-db", self.infiles.add_h_dict_infile.name,
                *shlex.split(self.parameters.reduce_options_parameter.value),
                self.infiles.receptor_infile.name,
            ],
            stdout_file_path=charged_receptor_full_h_file_path,
            ok_return_codes=(0, 1),  # 1: some flip optimizations were abandoned, but the output is complete
        )
        remove_reduce_annotations(charged_receptor_full_h_file_path)

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


def remove_reduce_annotations(pdb_file_path):
    """Remove the "new" flags reduce puts on added atoms and its USER records."""
    with open(pdb_file_path, "rb") as f:
        lines = f.read().split(b"\n")
    lines = [re.sub(rb"\s*new\s*", b"", line) for line in lines]
    lines = [line for line in lines if not line.startswith(b"USER")]
    with open(pdb_file_path, "wb") as f:
        f.write(b"\n".join(lines))
