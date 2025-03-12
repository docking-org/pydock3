import logging

from pydock3.blastermaster.util import BlasterStep, ProgramFilePaths
from pydock3.blastermaster.programs.visualization.create_VDW_DX import create_vdw_dx
from pydock3.blastermaster.programs.visualization.create_LigDeSolv_DX import create_ligdesolv_dx
from pydock3.blastermaster.programs.visualization.phi_to_dx import create_trim_electrostatics_dx

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class VisualizationStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        vdw_infile,
        vdw_bump_map_infile,
        lig_desolv_heavy_infile,
        trim_electrostatics_phi_infile,
        matching_spheres_infile,
        vdw_repulsive_dx_outfile,
        vdw_attractive_dx_outfile,
        vdw_dx_outfile,
        lig_desolv_dx_outfile,
        trim_electrostatics_dx_outfile,
        matching_spheres_outfile,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (vdw_infile, "vdw_infile", None),
                (vdw_bump_map_infile, "vdw_bump_map_infile", None),
                (lig_desolv_heavy_infile, "lig_desolv_heavy_infile", None),
                (trim_electrostatics_phi_infile, "trim_electrostatics_phi_infile", None),
                (matching_spheres_infile, "matching_spheres_infile", None)
            ],
            outfile_tuples=[
                (vdw_repulsive_dx_outfile, "vdw_repulsive_dx_outfile", None),
                (vdw_attractive_dx_outfile, "vdw_attractive_dx_outfile", None),
                (vdw_dx_outfile, "vdw_dx_outfile", None),
                (lig_desolv_dx_outfile, "lig_desolv_dx_outfile", None),
                (trim_electrostatics_dx_outfile, "trim_electrostatics_dx_outfile", None),
                (matching_spheres_outfile, "matching_spheres_outfile", None)
            ],
            parameter_tuples=[],
            program_file_path=None,
        )


    def vdw_visualization(self):
        create_vdw_dx(
            self.infiles.vdw_infile.path,
            self.infiles.vdw_bump_map_infile.path,
            self.outfiles.vdw_repulsive_dx_outfile.path,
            self.outfiles.vdw_attractive_dx_outfile.path,
            self.outfiles.vdw_dx_outfile.path,
        )

    def lig_desolv_visualization(self):
        create_ligdesolv_dx(
            self.infiles.lig_desolv_heavy_infile.path,
            self.outfiles.lig_desolv_dx_outfile.path
        )

    def trim_electrostatics_visualization(self):
        create_trim_electrostatics_dx(
            self.infiles.trim_electrostatics_phi_infile.path,
            self.outfiles.trim_electrostatics_dx_outfile.path
        )

    def visualize_matching_spheres(self):
        run_str = f"{ProgramFilePaths.DOSHOWSPH_PROGRAM_FILE_PATH} {self.infiles.matching_spheres_infile.name} 1 {self.outfiles.matching_spheres_outfile.name}"
        self.run_command(run_str)

    @BlasterStep.handle_run_func
    def run(self):
        self.vdw_visualization()
        self.lig_desolv_visualization()
        self.trim_electrostatics_visualization()
        self.visualize_matching_spheres()