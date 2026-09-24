import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.makebox import make_box
from pydock3.config import Parameter

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BoxGenerationStep(BlasterStep):

    # Default that can be overwritten in config
    MARGIN = Parameter("dock_files_generation.box_generation.margin", 10.0)

    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        box_outfile,
        margin_parameter=MARGIN,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
                (ligand_matching_spheres_infile, "ligand_matching_spheres_infile", None),
            ],
            outfile_tuples=[
                (box_outfile, "box_outfile", None),
            ],
            parameter_tuples=[
                (margin_parameter, "margin_parameter")
            ],
        )

    @BlasterStep.handle_run_func
    def run(self):
        """make box surrounding binding site"""
        make_box(
            self.infiles.ligand_matching_spheres_infile.path,
            self.infiles.charged_receptor_infile.path,
            self.outfiles.box_outfile.path,
            self.parameters.margin_parameter.value,
        )
