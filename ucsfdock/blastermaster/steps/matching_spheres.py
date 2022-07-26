#!/usr/bin/env python

# Ryan G. Coleman, Brian K. Shoichet Lab

import logging

from ucsfdock.blastermaster.util import ProgramFilePaths, BlasterStep
from ucsfdock.files import ProgramFile
from ucsfdock.blastermaster import pdb


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class MatchingSpheresGenerationStep(BlasterStep):

    MAX_NUM_SPHERES = 45
    DISTANCE_1 = 1.5
    DISTANCE_2 = 0.8

    def __init__(
        self,
        step_dir,
        charged_receptor_infile,
        ligand_matching_spheres_infile,
        all_spheres_infile,
        matching_spheres_outfile,
        covalent_use_parameter,
        covalent_residue_name_parameter,
        covalent_residue_num_parameter,
    ):
        super().__init__(step_dir=step_dir)

        #
        self.program_file = ProgramFile(
            path=ProgramFilePaths.MAKESPHERES3_PROGRAM_FILE_PATH
        )

        #
        self.process_infiles(
            (charged_receptor_infile, "charged_receptor_infile"),
            (ligand_matching_spheres_infile, "ligand_matching_spheres_infile"),
            (all_spheres_infile, "all_spheres_infile"),
        )

        #
        self.process_outfiles(
            (matching_spheres_outfile, "matching_spheres_outfile"),
        )

        #
        self.process_parameters(
            (covalent_use_parameter, "covalent_use_parameter"),
            (covalent_residue_name_parameter, "covalent_residue_name_parameter"),
            (covalent_residue_num_parameter, "covalent_residue_num_parameter"),
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run the makespheres3.cli.pl perl script to make low dielectric spheres"""

        # if this is a covalent run, export a matching spheres file based on the covalent residue
        if self.parameters.covalent_use_parameter.value:
            # output header coloring table for historical purposes
            with open(self.outfiles.matching_spheres_outfile.path, "w") as f:
                f.write("DOCK 5.2 ligand_atoms\n")
                f.write("positive                       (1)\n")
                f.write("negative                       (2)\n")
                f.write("acceptor                       (3)\n")
                f.write("donor                          (4)\n")
                f.write("ester_o                        (5)\n")
                f.write("amide_o                        (6)\n")
                f.write("neutral                        (7)\n")
                f.write("not_neutral                    (8)\n")
                f.write("positive_or_donor              (9)\n")
                f.write("negative_or_acceptor           (10)\n")
                f.write("neutral_or_acceptor_or_donor   (11)\n")
                f.write("donacc                         (12)\n")
                f.write("cluster     1   number of spheres in cluster   3\n")
                # read in receptor to get nucleophile coordinates
                pdb_h = pdb.PDBData(self.infiles.charged_receptor_infile.path, ignore_waters=False)

                if self.parameters.covalent_residue_name_parameter.value == "CYS":
                    p_1_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "CYS", "CA"
                        )
                    ]
                    p_2_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "CYS", "CB"
                        )
                    ]
                    p_3_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "CYS", "SG"
                        )
                    ]
                elif self.parameters.covalent_residue_name_parameter.value == "SER":
                    p_1_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "SER", "CA"
                        )
                    ]
                    p_2_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "SER", "CB"
                        )
                    ]
                    p_3_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "SER", "OG"
                        )
                    ]
                elif self.parameters.covalent_residue_name_parameter.value == "LYS":
                    p_1_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "LYS", "CD"
                        )
                    ]
                    p_2_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "LYS", "CE"
                        )
                    ]
                    p_3_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "LYS", "NZ"
                        )
                    ]
                elif self.parameters.covalent_residue_name_parameter.value == "TYR":
                    p_1_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "TYR", "CE1"
                        )
                    ]
                    p_2_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "TYR", "CZ"
                        )
                    ]
                    p_3_coords = pdb_h.coords[
                        pdb_h.get_index_by_residue_atom(
                            self.parameters.covalent_residue_num_parameter.value, "TYR", "OH"
                        )
                    ]
                else:
                    logger.exception(
                        f"Currently only supporting CYS, SER, LYS, TYR to prep other residues modify matching_spheres.py\ncovalent_residue_name given: {self.parameters.covalent_residue_name_parameter.value}"
                    )
                    raise

                f.write(
                    " 9001%10.5f%10.5f%10.5f   1.000    1 0  0\n"
                    % (p_1_coords[0], p_1_coords[1], p_1_coords[2])
                )
                f.write(
                    " 9002%10.5f%10.5f%10.5f   1.000    1 0  0\n"
                    % (p_2_coords[0], p_2_coords[1], p_2_coords[2])
                )
                f.write(
                    " 9003%10.5f%10.5f%10.5f   1.000    1 0  0\n"
                    % (p_3_coords[0], p_3_coords[1], p_3_coords[2])
                )

        else:
            # run
            run_str = f"{self.program_file.path} {self.DISTANCE_1} {self.DISTANCE_2} {self.MAX_NUM_SPHERES} {self.infiles.ligand_matching_spheres_infile.name} {self.infiles.all_spheres_infile.name} {self.infiles.charged_receptor_infile.name} {self.outfiles.matching_spheres_outfile.name}"
            self.run_command(run_str)
