import logging

from pydock3.blastermaster.util import ProgramFilePaths, BlasterStep
from pydock3.files import ProgramFile
from pydock3.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ChargedReceptorDeprotonationStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        charged_receptor_deprotonated_outfile,
        covalent_residue_num_parameter,
        covalent_residue_name_parameter,
        covalent_residue_atoms_parameter,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
            ],
            outfile_tuples=[
                (charged_receptor_deprotonated_outfile, "charged_receptor_deprotonated_outfile", None),
            ],
            parameter_tuples=[
                (covalent_residue_num_parameter, "covalent_residue_num_parameter"),
                (covalent_residue_name_parameter, "covalent_residue_name_parameter"),
                (covalent_residue_atoms_parameter, "covalent_residue_atoms_parameter"),
            ],
            program_file_path=ProgramFilePaths.REDUCE_PROGRAM_FILE_PATH,
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run REDUCE to produce a pdb with hydrogens.
        Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.
        then run script to remove nonpolar hydrogens & rename"""

        # deprotonate covalent residue
        logger.debug("Deprotonating covalent residue")
        pdb_h = pdb.PDBData(
            self.infiles.charged_receptor_infile.path,
            ignore_waters=False,
        )
        pdb_h.remove_protons_for_covalent_docking(
            self.parameters.covalent_residue_num_parameter.value,
            self.parameters.covalent_residue_name_parameter.value,
            self.parameters.covalent_residue_atoms_parameter.value,
        )
        pdb_h.write(self.outfiles.charged_receptor_deprotonated_outfile.path)
