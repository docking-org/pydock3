# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from pydock3.blastermaster.util import BlasterStep


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class ThinSpheresGenerationStep(BlasterStep):
    def __init__(
        self,
        step_dir,
        molecular_surface_infile,
        thin_spheres_outfile,
        distance_to_surface_parameter,
        penetration_parameter,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = None

        #
        self.process_infiles(
            (molecular_surface_infile, "molecular_surface_infile"),
        )

        #
        self.process_outfiles(
            (thin_spheres_outfile, "thin_spheres_outfile"),
        )

        #
        self.process_parameters(
            (distance_to_surface_parameter, "distance_to_surface_parameter"),
            (penetration_parameter, "penetration_parameter"),
        )

    @BlasterStep.handle_run_func
    def run(self):
        with open(self.infiles.molecular_surface_infile.path, "r") as f_in:
            with open(self.outfiles.thin_spheres_outfile.path, "w") as f_out:
                thin_spheres(f_in, f_out, distance=self.parameters.distance_to_surface_parameter.value, size=self.parameters.distance_to_surface_parameter.value + self.parameters.penetration_parameter.value)


def format_sphere_line(atom_num, sphere, size):
    """Format a line of a DOCK .sph file"""

    # FORTRAN FORMAT: (I5, 3F10.5, F8.3, I5, I2, I3)
    return "%5d%10.5f%10.5f%10.5f%8.3f%5d%2d%3d" % (atom_num, sphere[0], sphere[1], sphere[2], size, atom_num, 0, 0)


def thin_spheres(in_f, out_f, distance=1.8, size=1.9):
    """Generate delphi sphere pool near a given distance using rec.ms."""

    spheres = []
    for line in in_f:
        if line[40] == "S":
            splits = line.split()
            point = [float(line[13:21]), float(line[21:30]), float(line[30:39])]
            normal = [float(line[43:50]), float(line[50:57]), float(line[57:64]), float(line[64:71])]
            sphere = [p + distance * n for p, n in zip(point, normal)]

            # Handle odd chain id placement
            try:
                atom_num = int(splits[1])
            except ValueError:
                atom_num = int(splits[1][:-1])
            spheres.append((atom_num, sphere))

    out_f.write("cluster     0   number of spheres in cluster %5d\n" % len(spheres))
    for atom_num, sphere in spheres:
        out_f.write(format_sphere_line(atom_num, sphere, size=size) + "\n")
