# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorTransformationForLigandDesolvationNoThinSpheres(BlasterStep):
    def __init__(
        self,
        step_dir,
        charged_receptor_pdb_infile,
        charged_receptor_desolv_pdb_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (charged_receptor_pdb_infile, "charged_receptor_pdb_infile"),
        )

        #
        self.process_outfiles(
            (charged_receptor_desolv_pdb_outfile, "charged_receptor_desolv_pdb_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        #
        self.outfiles.charged_receptor_desolv_pdb_outfile.copy_from(self.infiles.charged_receptor_pdb_infile.path)


class ReceptorTransformationForLigandDesolvationYesThinSpheres(BlasterStep):
    def __init__(
        self,
        step_dir,
        charged_receptor_pdb_infile,
        close_spheres_desolv_pdb_infile,
        charged_receptor_desolv_pdb_outfile,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (charged_receptor_pdb_infile, "charged_receptor_pdb_infile"),
            (close_spheres_desolv_pdb_infile, "close_spheres_desolv_pdb_infile"),
        )

        #
        self.process_outfiles(
            (charged_receptor_desolv_pdb_outfile, "charged_receptor_desolv_pdb_outfile"),
        )

        #
        self.process_parameters()

    @BlasterStep.handle_run_func
    def run(self):
        #
        self.outfiles.charged_receptor_desolv_pdb_outfile.copy_from(self.infiles.charged_receptor_pdb_infile.path)

        #
        run_str = f"cat {self.infiles.close_spheres_desolv_pdb_infile.path} | sed -e 's/ C   SPH/ X   SPH/g' >> {self.outfiles.charged_receptor_desolv_pdb_outfile.path}"
        self.run_command(run_str)
