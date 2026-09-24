import os
import re
import logging

from pydock3.blastermaster.util import program_path, BlasterStep
from pydock3.files import File
from pydock3.blastermaster import phi
from pydock3.config import Parameter


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_scale_and_center(log_file):
    """reads a logfile, finds the scale and center lines, returns them"""
    with open(log_file.path, "r") as f:
        lines = ["sizing=scale"]
        for line in f:
            for match_value in ["scale=", "xcen=", "ycen=", "zcen="]:
                try:
                    if match_value in line:
                        lines.append(line.strip())
                except ValueError:
                    pass
        return lines


class ElectrostaticsGridGenerationStepNoThinSpheres(BlasterStep):

    QNIFFT_PARAMETERS_FILE_NAME: str = "qnifft.parm"

    # Defaults that can be overwritten in the config
    GRID_SIZE = Parameter("dock_files_generation.electrostatics_grid_gen.grid_size", 193)
    USE_RECEPTOR_BOX = Parameter("dock_files_generation.electrostatics_grid_gen.use_receptor_box", False)

    def __init__(
        self,
        working_dir,
        receptor_low_dielectric_pdb_infile,
        charge_infile,
        radius_infile,
        delphi_infile,
        box_infile,
        electrostatics_phi_outfile,
        electrostatics_pdb_outfile,
        electrostatics_trim_phi_outfile,
        electrostatics_phi_size_outfile,
        grid_size_parameter=GRID_SIZE,
        use_receptor_box_parameter=USE_RECEPTOR_BOX,
        extra_parameters=None,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_low_dielectric_pdb_infile, "receptor_low_dielectric_pdb_infile", None),
                (charge_infile, "charge_infile", None),
                (radius_infile, "radius_infile", None),
                (delphi_infile, "delphi_infile", None),
                (box_infile, "box_infile", None),
            ],
            outfile_tuples=[
                (electrostatics_phi_outfile, "electrostatics_phi_outfile", None),
                (electrostatics_pdb_outfile, "electrostatics_pdb_outfile", None),
                (electrostatics_trim_phi_outfile, "electrostatics_trim_phi_outfile", None),
                (electrostatics_phi_size_outfile, "electrostatics_phi_size_outfile", None),
            ],
            parameter_tuples=[
                (grid_size_parameter, "grid_size_parameter"),
                (use_receptor_box_parameter, "use_receptor_box_parameter")
            ],
        )

        # misc.
        self.extra_parameters = extra_parameters

    @BlasterStep.handle_run_func
    def run(self):
        """run the electrostatics program qnifft"""

        if self.extra_parameters is None:
            self.extra_parameters = []

        # first make parm file
        qnifft_parameters_file = File(
            path=os.path.join(self.step_dir.path, self.QNIFFT_PARAMETERS_FILE_NAME)
        )
        with open(qnifft_parameters_file.path, "w", newline="\n") as f:
            with open(self.infiles.delphi_infile.path, "r") as g:
                for line in g:
                    f.write(line)
            f.write(f"grid={str(self.parameters.grid_size_parameter.value)}\n")
            f.write(f"charge={self.infiles.charge_infile.name}\n")
            f.write(f"radius={self.infiles.radius_infile.name}\n")
            f.write(
                f"pdb_input={self.infiles.receptor_low_dielectric_pdb_infile.name}\n"
            )
            f.write(
                f"pdb_output_file={self.outfiles.electrostatics_pdb_outfile.name}\n"
            )
            f.write(
                f"phi_output_file={self.outfiles.electrostatics_phi_outfile.name}\n"
            )
            if (
                self.parameters.use_receptor_box_parameter.value
            ):  # means we are running the whole protein + 8 angstroms
                f.write("border=15\n")
            if self.extra_parameters:
                for extra_param in self.extra_parameters:
                    f.write(f"{extra_param}\n")
            else:
                f.write("sizing=border border_solvent=10.\n")

        #
        self.log_parameters_file(qnifft_parameters_file)

        # run program
        self.run_program([program_path("qnifft"), qnifft_parameters_file.name])

        #
        phi_size, new_center = phi.trim(
            input_phi_file=self.outfiles.electrostatics_phi_outfile,
            box_file=self.infiles.box_infile,
            output_phi_file=self.outfiles.electrostatics_trim_phi_outfile,
        )

        #
        with open(self.outfiles.electrostatics_phi_size_outfile.path, "w") as f:
            f.write(str(phi_size))


class ElectrostaticsGridGenerationStepYesThinSpheres(BlasterStep):

    QNIFFT_PARAMETERS_FILE_NAME: str = "qnifft.parm"

    # Defaults that can be overwritten in the config
    GRID_SIZE = Parameter("dock_files_generation.electrostatics_grid_gen.grid_size", 193)
    USE_RECEPTOR_BOX = Parameter("dock_files_generation.electrostatics_grid_gen.use_receptor_box", False)

    def __init__(
        self,
        working_dir,
        receptor_low_dielectric_pdb_infile,
        charge_infile,
        radius_infile,
        delphi_infile,
        box_infile,
        electrostatics_phi_outfile,
        electrostatics_pdb_outfile,
        electrostatics_trim_phi_outfile,
        electrostatics_phi_size_outfile,
        thin_spheres_elec_distance_to_surface_parameter,
        thin_spheres_elec_penetration_parameter,
        grid_size_parameter=GRID_SIZE,
        use_receptor_box_parameter=USE_RECEPTOR_BOX,
        extra_parameters=None,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (receptor_low_dielectric_pdb_infile, "receptor_low_dielectric_pdb_infile", None),
                (charge_infile, "charge_infile", None),
                (radius_infile, "radius_infile", None),
                (delphi_infile, "delphi_infile", None),
                (box_infile, "box_infile", None),
            ],
            outfile_tuples=[
                (electrostatics_phi_outfile, "electrostatics_phi_outfile", None),
                (electrostatics_pdb_outfile, "electrostatics_pdb_outfile", None),
                (electrostatics_trim_phi_outfile, "electrostatics_trim_phi_outfile", None),
                (electrostatics_phi_size_outfile, "electrostatics_phi_size_outfile", None),
            ],
            parameter_tuples=[
                (thin_spheres_elec_distance_to_surface_parameter, "thin_spheres_elec_distance_to_surface_parameter"),
                (thin_spheres_elec_penetration_parameter, "thin_spheres_elec_penetration_parameter"),
                (grid_size_parameter, "grid_size_parameter"),
                (use_receptor_box_parameter, "use_receptor_box_parameter")
            ],
        )

        # misc.
        self.extra_parameters = extra_parameters

    @BlasterStep.handle_run_func
    def run(self):
        """run the electrostatics program qnifft"""

        if self.extra_parameters is None:
            self.extra_parameters = []

        # first make parm file
        qnifft_parameters_file = File(
            path=os.path.join(self.step_dir.path, self.QNIFFT_PARAMETERS_FILE_NAME)
        )
        with open(qnifft_parameters_file.path, "w", newline="\n") as f:
            with open(self.infiles.delphi_infile.path, "r") as g:
                for line in g:
                    f.write(line)
            f.write(f"grid={str(self.parameters.grid_size_parameter.value)}\n")
            f.write(f"charge={self.infiles.charge_infile.name}\n")
            f.write(f"radius={self.infiles.radius_infile.name}\n")
            f.write(
                f"pdb_input={self.infiles.receptor_low_dielectric_pdb_infile.name}\n"
            )
            f.write(
                f"pdb_output_file={self.outfiles.electrostatics_pdb_outfile.name}\n"
            )
            f.write(
                f"phi_output_file={self.outfiles.electrostatics_phi_outfile.name}\n"
            )
            if (
                self.parameters.use_receptor_box_parameter.value
            ):  # means we are running the whole protein + 8 angstroms
                f.write("border=15\n")
            if self.extra_parameters:
                for extra_param in self.extra_parameters:
                    f.write(f"{extra_param}\n")
            else:
                f.write("sizing=border border_solvent=10.\n")

        #
        self.log_parameters_file(qnifft_parameters_file)

        # the thin spheres get the radius given by their distance to the surface
        sphere_radius = (
            self.parameters.thin_spheres_elec_distance_to_surface_parameter.value
            + self.parameters.thin_spheres_elec_penetration_parameter.value
        )
        with open(self.infiles.radius_infile.path, "rb") as f:
            radii = f.read()
        with open(self.infiles.radius_infile.path, "wb") as f:
            f.write(re.sub(rb"c     sph   1.90", b"c     sph   %3.2f" % sphere_radius, radii))

        # run program
        self.run_program([program_path("qnifft"), qnifft_parameters_file.name])

        #
        phi_size, new_center = phi.trim(
            input_phi_file=self.outfiles.electrostatics_phi_outfile,
            box_file=self.infiles.box_infile,
            output_phi_file=self.outfiles.electrostatics_trim_phi_outfile,
        )

        #
        with open(self.outfiles.electrostatics_phi_size_outfile.path, "w") as f:
            f.write(str(phi_size))
