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
            program_file_path=None,
        )

    @BlasterStep.handle_run_func
    def run(self):
        run_str = f"cat {self.infiles.charged_receptor_infile.name} {self.infiles.spheres_pdb_infile.name} > {self.outfiles.receptor_low_dielectric_pdb_outfile.name}"
        self.run_command(run_str)
