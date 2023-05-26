import os
import tarfile
from typing import List


class RetrospectiveDataset(object):
    SUPPORTED_EXTENSIONS = ['db2', 'db2.gz']

    def __init__(self, positives_tgz_file_path: str, negatives_tgz_file_path: str, positives_dir_path: str, negatives_dir_path: str):
        #
        if not os.path.isfile(positives_tgz_file_path):
            raise Exception(f"Tarball `{positives_tgz_file_path}` does not exist.")
        if not os.path.isfile(negatives_tgz_file_path):
            raise Exception(f"Tarball `{negatives_tgz_file_path}` does not exist.")

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
        self.num_positives = self._extract_tarball_and_count_files(self.positives_tgz_file_path, self.SUPPORTED_EXTENSIONS, self.positives_dir_path)
        self.num_negatives = self._extract_tarball_and_count_files(self.negatives_tgz_file_path, self.SUPPORTED_EXTENSIONS, self.negatives_dir_path)

        #
        if self.num_positives == 0:
            raise Exception(f"No positives found in tarball `{self.positives_tgz_file_path}`. Expected files with extensions: `{self.SUPPORTED_EXTENSIONS}`")
        if self.num_negatives == 0:
            raise Exception(f"No negatives found in tarball `{self.negatives_tgz_file_path}`. Expected files with extensions: `{self.SUPPORTED_EXTENSIONS}`")

    @staticmethod
    def _list_tarball_contents(tarball):
        with tarfile.open(tarball, "r:gz") as tar:
            return sorted([tarinfo.name for tarinfo in tar])

    @staticmethod
    def _list_directory_contents(directory):
        return sorted([f for root, dirs, files in os.walk(directory) for f in files])

    @staticmethod
    def _check_that_directory_and_tarball_match(directory, tarball):
        # Get contents
        tarball_contents = RetrospectiveDataset._list_tarball_contents(tarball)
        directory_contents = RetrospectiveDataset._list_directory_contents(directory)

        # Remove extensions
        tarball_contents = [f[:-3] if f.endswith(".gz") else f for f in tarball_contents]
        directory_contents = [f[:-3] if f.endswith(".gz") else f for f in directory_contents]

        return tarball_contents == directory_contents

    @staticmethod
    def _extract_tarball_and_count_files(tarball_path: str, supported_extensions: List[str], extract_dir_path: str):
        if not os.path.isfile(tarball_path):
            raise Exception(f"Tarball `{tarball_path}` does not exist.")

        if not os.path.isdir(extract_dir_path):
            raise Exception(f"Directory `{extract_dir_path}` does not exist.")

        file_count = 0
        with tarfile.open(tarball_path, 'r:gz') as tar:
            #
            for member in tar.getmembers():
                if member.isfile() and any([member.name.lower().endswith(f".{ext}") for ext in supported_extensions]):
                    file_count += 1
                elif member.isfile():
                    raise Exception(
                        f"File `{member.name}` in tarball `{tarball_path}` has an unsupported extension. Extension must be one of: `{supported_extensions}`")

            #
            if len(os.listdir(extract_dir_path)) > 0:  # dir has already been extracted
                if not RetrospectiveDataset._check_that_directory_and_tarball_match(extract_dir_path, tarball_path):
                    raise Exception(f"Contents of extract directory `{extract_dir_path}` do not match those of tarball `{tarball_path}`. If you made a change to the tarball, you must delete the extract directory and re-run.")
            else:  # dir has not been extracted yet, so extract it
                tar.extractall(path=extract_dir_path)

        if file_count == 0:
            raise Exception(f"No files with supported extensions `{supported_extensions}` found in tarball `{tarball_path}`.")

        return file_count
