import logging
from copy import deepcopy
import random

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
            ],
            program_file_path=None,
        )

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
        seed = get_hexdigest_of_persistent_md5_hash_of_tuple(tuple(sphere_hashes + [self.outfiles.perturbed_matching_spheres_outfile.name]))
        random.seed(seed)

        # perturb all spheres in file
        new_spheres = []
        for sphere in spheres:
            new_sphere = deepcopy(sphere)
            max_deviation = float(self.parameters.max_deviation_angstroms.value)
            perturbation_xyz = tuple(
                [
                    random.uniform(
                        -max_deviation,
                        max_deviation,
                    )
                    for _ in range(3)
                ]
            )
            new_sphere.X += perturbation_xyz[0]
            new_sphere.Y += perturbation_xyz[1]
            new_sphere.Z += perturbation_xyz[2]
            new_spheres.append(new_sphere)

        # write perturbed spheres to new matching spheres file
        write_sph(self.outfiles.perturbed_matching_spheres_outfile.path, new_spheres)
