import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.sphere_files import sph_to_pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class SpheresToPDBConversionStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        sph_infile,
        pdb_outfile,
    ):
        #
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (sph_infile, "sph_infile", None),
            ],
            outfile_tuples=[
                (pdb_outfile, "pdb_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        """convert spheres to pdb file"""

        sph_to_pdb(self.infiles.sph_infile.path, 1, self.outfiles.pdb_outfile.path)
