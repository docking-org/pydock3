#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorTransformationForElectrostatics(BlasterStep):
    def __init__(
        self,
        step_dir,
        charged_receptor_infile,
        spheres_pdb_infile,
        receptor_low_dielectric_pdb_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (charged_receptor_infile, "charged_receptor_infile"),
            (spheres_pdb_infile, "spheres_pdb_infile"),
        )

        #
        self.process_outfiles(
            (receptor_low_dielectric_pdb_outfile, "receptor_low_dielectric_pdb_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        run_str = f"cat {self.infiles.charged_receptor_infile.name} {self.infiles.spheres_pdb_infile.name} > {self.outfiles.receptor_low_dielectric_pdb_outfile.name}"
        self.run_command(run_str)
