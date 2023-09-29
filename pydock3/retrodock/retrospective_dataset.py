import os
import tarfile
from typing import List, Generator

from pydock3.files import TarballFile, DB2File


class RetrospectiveDataset(object):
    SUPPORTED_EXTENSIONS = ['db2', 'db2.gz']

    def __init__(self, positives_tgz_file_path: str, negatives_tgz_file_path: str, positives_dir_path: str, negatives_dir_path: str):
        #
        if not os.path.isfile(positives_tgz_file_path):
            raise Exception(f"Tarball `{positives_tgz_file_path}` does not exist.")
        if not os.path.isfile(negatives_tgz_file_path):
            raise Exception(f"Tarball `{negatives_tgz_file_path}` does not exist.")

        #
        self._validate_tarball_files(positives_tgz_file_path)
        self._validate_tarball_files(negatives_tgz_file_path)

        # Create directories if they don't exist
        for dir_path in [positives_dir_path, negatives_dir_path]:
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)

        #
        self.positives_tgz_file_path = positives_tgz_file_path
        self.negatives_tgz_file_path = negatives_tgz_file_path

        #
        self.positives_dir_path = positives_dir_path
        self.negatives_dir_path = negatives_dir_path

        #
        self.num_db2_files_in_positive_class = len(list(TarballFile(self.positives_tgz_file_path).iterate_over_tarball_member_files()))
        self.num_db2_files_in_negative_class = len(list(TarballFile(self.negatives_tgz_file_path).iterate_over_tarball_member_files()))

        #
        if self.num_db2_files_in_positive_class == 0:
            raise Exception(f"No positives found in tarball `{self.positives_tgz_file_path}`. Expected files with extensions: `{self.SUPPORTED_EXTENSIONS}`")
        if self.num_db2_files_in_negative_class == 0:
            raise Exception(f"No negatives found in tarball `{self.negatives_tgz_file_path}`. Expected files with extensions: `{self.SUPPORTED_EXTENSIONS}`")

        #
        self._extract_tarball(self.positives_tgz_file_path, self.positives_dir_path)
        self._extract_tarball(self.negatives_tgz_file_path, self.negatives_dir_path)

        #
        self.num_molecules_in_positive_class = len(list(set([DB2File(os.path.join(self.positives_dir_path, file.name.lstrip('./'))).get_molecule_name() for file in TarballFile(self.positives_tgz_file_path).iterate_over_tarball_member_files()])))
        self.num_molecules_in_negative_class = len(list(set([DB2File(os.path.join(self.negatives_dir_path, file.name.lstrip('./'))).get_molecule_name() for file in TarballFile(self.negatives_tgz_file_path).iterate_over_tarball_member_files()])))

    def _validate_tarball_files(self, tarball_path: str) -> None:
        file_count = 0
        for file in TarballFile(tarball_path).iterate_over_tarball_member_files():
            if (any([file.name.lower().endswith(f".{ext}") for ext in self.SUPPORTED_EXTENSIONS])):
                file_count += 1
            else:
                raise Exception(f"File `{file.name}` in tarball `{tarball_path}` has an unsupported extension. Extension must be one of: `{self.SUPPORTED_EXTENSIONS}`")

        if file_count == 0:
            raise Exception(f"No files with supported extensions `{self.SUPPORTED_EXTENSIONS}` found in tarball `{tarball_path}`.")

    @staticmethod
    def _extract_tarball(tarball_path: str, extraction_dir_path: str):
        def _check_that_extraction_directory_and_tarball_match(extraction_dir_path: str, tarball_path: str):
            # Get the set of file names (including directory structure) in the tarball
            tarball_files = set(
                file.name.lstrip('./') for file in TarballFile(tarball_path).iterate_over_tarball_member_files())

            # Recursively traverse through all subdirectories and files in the extraction directory
            # and calculate the relative path of each file with respect to the extraction directory
            extracted_files = set(os.path.relpath(
                os.path.join(root, file), extraction_dir_path)
                for root, dirs, files in os.walk(extraction_dir_path)
                for file in files
            )

            return tarball_files == extracted_files

        if len(os.listdir(extraction_dir_path)) > 0:  # dir has already been extracted
            if not _check_that_extraction_directory_and_tarball_match(extraction_dir_path, tarball_path):
                raise Exception(f"Contents of extract directory `{extraction_dir_path}` do not match those of tarball `{tarball_path}`. If you made a change to the tarball, you must delete the extract directory and re-run.")
        else:  # dir has not been extracted yet, so extract it
            TarballFile(tarball_path).extract(extraction_dir_path)
