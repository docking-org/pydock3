import filecmp
import os

from pydock3.dockopt.dockopt import Dockopt


def test_dockopt(tmp_path):
    #
    dockopt = Dockopt()

    #
    tmp_job_path = os.path.join(tmp_path, "dockopt_job")
    dockopt.init(job_dir_path=tmp_job_path)

    # replace dockopt_config.yaml, input files, and retrospective docking data with test versions
    # TODO

    #
    dockopt.run(
        scheduler=None,
        job_dir_path=tmp_job_path,
    )

    # validate output files exist
    # TODO
    for output_file in output_files:
        assert output_file.exists

    # compare output files with test versions
    # TODO
    for output_file in output_files:
        assert filecmp.cmp(output_file.path, test_output_file.path)