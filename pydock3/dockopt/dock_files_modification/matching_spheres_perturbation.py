import logging
from copy import deepcopy
import random
import numpy as np

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.thinspheres.sph_lib import read_sph, write_sph
from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class MatchingSpheresPerturbationStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        matching_spheres_infile,
        perturbed_matching_spheres_outfile,
        max_deviation_angstroms_parameter,
        perturb_xtal_spheres_parameter,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (matching_spheres_infile, "matching_spheres_infile", None),
            ],
            outfile_tuples=[
                (perturbed_matching_spheres_outfile, "perturbed_matching_spheres_outfile", None),
            ],
            parameter_tuples=[
                (max_deviation_angstroms_parameter, "max_deviation_angstroms_parameter"),
                (perturb_xtal_spheres_parameter, "perturb_xtal_spheres_parameter")
            ],
            program_file_path=None,
        )

    # UNIFORMLY samples a distance (0,max_r) to move the center
    # Not the previous algorithm that moved a sphere inside a cube ;) 
    def move_sphere_center(self,x,y,z,max_r):
        # Sample distance uniformly in [0,max_r]
        d = np.random.uniform(0, max_r)

        # Sample a random direction on the sphere
        phi = np.random.uniform(0, 2 * np.pi) # Azimuthal angle
        cos_theta = np.random.uniform(-1, 1) # Cosine of polar angle
        theta = np.arccos(cos_theta)
        
        # Convert to cartesian
        dx = d * np.sin(theta) * np.cos(phi)
        dy = d * np.sin(theta) * np.sin(phi)
        dz = d * np.cos(theta)

        new_x = x + dx
        new_y = y + dy
        new_z = z + dz

        return new_x,new_y,new_z

    @BlasterStep.handle_run_func
    def run(self):
        """#TODO"""

        #
        spheres = read_sph(
            self.infiles.matching_spheres_infile.path,
            chosen_cluster="A",
            color="A",
        )

        # set random seed based on spheres and outfile name for reproducibility
        sphere_hashes = [get_hexdigest_of_persistent_md5_hash_of_tuple((sphere.index, sphere.X, sphere.Y, sphere.Z, sphere.radius, sphere.atomnum, sphere.critical_cluster, sphere.sphere_color)) for sphere in spheres]
        seed = get_hexdigest_of_persistent_md5_hash_of_tuple(tuple(sphere_hashes + [self.outfiles.perturbed_matching_spheres_outfile.name, self.parameters.max_deviation_angstroms_parameter.value]))
        seed = int(seed, 16)  # Convert hex string to an integer
        np.random.seed(seed % (2**32))  # Keep within valid range

        # perturb spheres
        new_spheres = []
        for sphere in spheres:
            # if the sphere is from xtal-lig then the radius is 0.5. Then check if we want to perturb it
            if sphere.radius == 0.5 and not self.parameters.perturb_xtal_spheres_parameter.value:
                new_spheres.append(deepcopy(sphere))
                continue
            new_sphere = deepcopy(sphere)
            max_deviation = float(self.parameters.max_deviation_angstroms_parameter.value)
            new_sphere.X, new_sphere.Y, new_sphere.Z = self.move_sphere_center(new_sphere.X, new_sphere.Y, new_sphere.Z, max_deviation)
            new_spheres.append(new_sphere)

        # write perturbed spheres to new matching spheres file
        write_sph(self.outfiles.perturbed_matching_spheres_outfile.path, new_spheres)
