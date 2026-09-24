import logging

from pydock3.blastermaster.util import program_path, BlasterStep
from pydock3.files import File
from pydock3.config import Parameter 


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class MolecularSurfaceGenerationStep(BlasterStep):

    # Default value that can be overwritten in config
    DENSITY = Parameter("dock_files_generation.molecular_surface_density", 5.0)

    class MandatoryFileNames:
        RADII_FILE_NAME = "radii"

    def __init__(
        self,
        working_dir,
        charged_receptor_infile,
        binding_site_residues_infile,
        radii_infile,
        molecular_surface_outfile,
        molecular_surface_density_parameter=DENSITY,
    ):
        super().__init__(
            working_dir=working_dir,
            infile_tuples=[
                (charged_receptor_infile, "charged_receptor_infile", None),
                (binding_site_residues_infile, "binding_site_residues_infile", None),
                (radii_infile, "radii_infile", self.MandatoryFileNames.RADII_FILE_NAME), # dms reads the elements and radii from a file in the current directory called 'radii'
            ],
            outfile_tuples=[
                (molecular_surface_outfile, "molecular_surface_outfile", None),
            ],
            parameter_tuples=[
                (molecular_surface_density_parameter, "molecular_surface_density_parameter")
            ],
        )

    @BlasterStep.handle_run_func
    def run(self):
        """run the dms program to produce molecular surface points
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html
        """
        # if you have waters, DMS crashes. so, just take them out of any file DMS
        # reads. this should allow statically placed waters to work.
        # thanks to joel karpiak for finding these errors!

        # removing waters from receptors
        charged_receptor_no_waters_file = File(
            path=f"{self.infiles.charged_receptor_infile.path}.dms"
        )
        remove_lines_containing(b"HOH", self.infiles.charged_receptor_infile.path, charged_receptor_no_waters_file.path)

        # removing waters from binding site
        binding_site_residues_no_waters_file = File(
            path=f"{self.infiles.binding_site_residues_infile.path}.dms"
        )
        remove_lines_containing(b"HOH", self.infiles.binding_site_residues_infile.path, binding_site_residues_no_waters_file.path)

        #
        self.run_program([
            program_path("dms"), charged_receptor_no_waters_file.name,
            "-a", "-d", self.parameters.molecular_surface_density_parameter.value,
            "-i", binding_site_residues_no_waters_file.name,
            "-g", self.log_file.name, "-p", "-n",
            "-o", self.outfiles.molecular_surface_outfile.name,
        ])


def remove_lines_containing(text, in_file_path, out_file_path):
    """Copy a file without the lines containing `text` (bytes), like `grep -a -v`."""
    with open(in_file_path, "rb") as f:
        lines = f.read().split(b"\n")
    if lines[-1] == b"":  # the file ended with a newline
        lines.pop()
    with open(out_file_path, "wb") as f:
        f.writelines(line + b"\n" for line in lines if text not in line)
