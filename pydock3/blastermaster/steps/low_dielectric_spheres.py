import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.makespheres import make_low_dielectric_spheres
from pydock3.config import Parameter

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LowDielectricSpheresSelectionStep(BlasterStep):

    # Default that can be overwritten in the config
    MIN_NUM_SPHERES = Parameter("dock_files_generation.low_dielectric_sphere_selection.min_num_spheres", 25)

    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        all_spheres_infile,
        low_dielectric_spheres_outfile,
        min_num_spheres_parameter=MIN_NUM_SPHERES
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
                (ligand_matching_spheres_infile, "ligand_matching_spheres_infile", None),
                (all_spheres_infile, "all_spheres_infile", None),
            ],
            outfile_tuples=[
                (low_dielectric_spheres_outfile, "low_dielectric_spheres_outfile", None),
            ],
            parameter_tuples=[
                (min_num_spheres_parameter, "min_num_spheres_parameter")
            ],
        )


    @BlasterStep.handle_run_func
    def run(self):
        """make low dielectric spheres"""
        make_low_dielectric_spheres(
            self.infiles.ligand_matching_spheres_infile.path,
            self.infiles.all_spheres_infile.path,
            self.infiles.charged_receptor_infile.path,
            self.outfiles.low_dielectric_spheres_outfile.path,
            self.parameters.min_num_spheres_parameter.value,
        )
