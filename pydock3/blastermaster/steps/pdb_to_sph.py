import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.sphere_files import pdb_to_sph


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LigandPDBToSpheresConversionStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        pdb_infile,
        sph_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (pdb_infile, "pdb_infile", None),
            ],
            outfile_tuples=[
                (sph_outfile, "sph_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        """turn the ligand into spheres"""
        pdb_to_sph(self.infiles.pdb_infile.path, self.outfiles.sph_outfile.path)
