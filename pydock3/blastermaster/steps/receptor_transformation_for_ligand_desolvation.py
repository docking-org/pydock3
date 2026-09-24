import logging

from pydock3.blastermaster.util import BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ReceptorTransformationForLigandDesolvationNoThinSpheres(BlasterStep):
    def __init__(
        self,
        working_dir,
        charged_receptor_pdb_infile,
        charged_receptor_desolv_pdb_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_pdb_infile, "charged_receptor_pdb_infile", None),
            ],
            outfile_tuples=[
                (charged_receptor_desolv_pdb_outfile, "charged_receptor_desolv_pdb_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        #
        self.outfiles.charged_receptor_desolv_pdb_outfile.copy_from(
            self.infiles.charged_receptor_pdb_infile.path
        )


class ReceptorTransformationForLigandDesolvationYesThinSpheres(BlasterStep):
    def __init__(
        self,
        working_dir,
        charged_receptor_pdb_infile,
        close_spheres_desolv_pdb_infile,
        charged_receptor_desolv_pdb_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_pdb_infile, "charged_receptor_pdb_infile", None),
                (close_spheres_desolv_pdb_infile, "close_spheres_desolv_pdb_infile", None),
            ],
            outfile_tuples=[
                (charged_receptor_desolv_pdb_outfile, "charged_receptor_desolv_pdb_outfile", None),
            ],
            parameter_tuples=[],
        )

    @BlasterStep.handle_run_func
    def run(self):
        #
        self.outfiles.charged_receptor_desolv_pdb_outfile.copy_from(
            self.infiles.charged_receptor_pdb_infile.path
        )

        # append the spheres, as atom type X
        with open(self.infiles.close_spheres_desolv_pdb_infile.path, "rb") as f:
            spheres = f.read().replace(b" C   SPH", b" X   SPH")
        with open(self.outfiles.charged_receptor_desolv_pdb_outfile.path, "ab") as f:
            f.write(spheres)
