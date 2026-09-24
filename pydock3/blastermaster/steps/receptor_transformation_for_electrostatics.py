import logging

from pydock3.blastermaster.util import BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorTransformationForElectrostatics(BlasterStep):
    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        spheres_pdb_infile,
        receptor_low_dielectric_pdb_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
                (spheres_pdb_infile, "spheres_pdb_infile", None),
            ],
            outfile_tuples=[
                (receptor_low_dielectric_pdb_outfile, "receptor_low_dielectric_pdb_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        with open(self.outfiles.receptor_low_dielectric_pdb_outfile.path, "wb") as f_out:
            for infile in (self.infiles.charged_receptor_infile, self.infiles.spheres_pdb_infile):
                with open(infile.path, "rb") as f_in:
                    f_out.write(f_in.read())
