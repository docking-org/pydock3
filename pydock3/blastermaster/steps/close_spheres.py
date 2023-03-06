import logging

from pydock3.blastermaster.util import BlasterStep
from pydock3.blastermaster.programs.thinspheres import sph_lib, pdb_lib

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class CloseSpheresGenerationStep(BlasterStep):
    def __init__(
        self,
        working_dir,
        ligand_infile,
        thin_spheres_infile,
        close_spheres_outfile,
        distance_to_surface_parameter,
        penetration_parameter,
        distance_to_ligand_parameter,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (ligand_infile, "ligand_infile", None),
                (thin_spheres_infile, "thin_spheres_infile", None),
            ],
            outfile_tuples=[
                (close_spheres_outfile, "close_spheres_outfile", None),
            ],
            parameter_tuples=[
                (distance_to_surface_parameter, "distance_to_surface_parameter"),
                (penetration_parameter, "penetration_parameter"),
                (distance_to_ligand_parameter, "distance_to_ligand_parameter"),
            ],
            program_file_path=None,
        )

    @BlasterStep.handle_run_func
    def run(self):
        spheres_list = sph_lib.read_sph(self.infiles.thin_spheres_infile.path, "A", "A")
        pdb_list = pdb_lib.read_pdb(self.infiles.ligand_infile.path)
        spheres_list = distance_sph_pdb(
            spheres_list, pdb_list, self.parameters.distance_to_ligand_parameter.value
        )
        radius = (
            self.parameters.distance_to_surface_parameter.value
            + self.parameters.penetration_parameter.value
        )
        spheres_list = trim_sph(spheres_list, radius)
        sph_lib.write_sph(self.outfiles.close_spheres_outfile.path, spheres_list)


def trim_sph(sph_list, sph_rad):
    for i in range(len(sph_list) - 1):
        if sph_list[i][1]:
            for j in range(i + 1, len(sph_list)):
                if sph_list[j][1]:
                    dist = (
                        (sph_list[i][0].X - sph_list[j][0].X) ** 2
                        + (sph_list[i][0].Y - sph_list[j][0].Y) ** 2
                        + (sph_list[i][0].Z - sph_list[j][0].Z) ** 2
                    )
                    if dist <= (sph_rad**2.0) / 2.0:
                        sph_list[j][1] = False

    final_sph_list = []
    for sph in sph_list:
        if sph[1]:
            final_sph_list.append(sph[0])

    return final_sph_list


def distance_sph_pdb(spheres, pdb_atoms, distance):
    sph_list = []
    for sph in spheres:
        for atom in pdb_atoms:
            d2 = (atom.X - sph.X) ** 2 + (atom.Y - sph.Y) ** 2 + (atom.Z - sph.Z) ** 2
            if d2 < float(distance) ** 2.0:
                sph_list.append([sph, True])
                break

    return sph_list
