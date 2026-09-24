import os
import logging

from pydock3.blastermaster.util import program_path, BlasterStep
from pydock3.files import File


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BindingSiteSpheresGenerationStep(BlasterStep):

    INSPH_FILE_NAME = "INSPH"  # Input file needs to have this name for sphgen to run
    LINES_2_THROUGH_6 = (
        "R\nX\n0.\n5.0\n1.4"  # line1 is input file (rec.ms), line7 is output file (sph)
    )

    def __init__(
        self,
        working_dir,
        molecular_surface_infile,
        spheres_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (molecular_surface_infile, "molecular_surface_infile", None),
            ],
            outfile_tuples=[
                (spheres_outfile, "spheres_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run the sphgen program to produce initial set of spheres"""

        # make input file
        sphgen_input_file = File(
            path=os.path.join(self.step_dir.path, self.INSPH_FILE_NAME)
        )
        with open(sphgen_input_file.path, "w", newline="\n") as f:
            f.write(f"{self.infiles.molecular_surface_infile.name}\n")  # input file
            f.write(f"{self.LINES_2_THROUGH_6}\n")  # see above
            f.write(f"{self.outfiles.spheres_outfile.name}\n")  # output file

        # run
        self.run_program([program_path("sphgen")])

        # remove the first line of the output
        with open(self.outfiles.spheres_outfile.path, "rb") as f:
            _, _, rest = f.read().partition(b"\n")
        with open(self.outfiles.spheres_outfile.path, "wb") as f:
            f.write(rest)
