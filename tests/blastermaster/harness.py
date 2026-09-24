"""Run blastermaster jobs on the test fixtures (shared by the tests and `make_controls.py`)."""
import os
import shutil
from contextlib import contextmanager
from pathlib import Path

import oyaml as yaml

from pydock3.blastermaster.blastermaster import Blastermaster, load_steps
from pydock3.files import INDOCK_FILE_NAME

HERE = Path(__file__).resolve().parent
FIXTURE_FILES = [HERE.parent / "rec.pdb", HERE.parent / "xtal-lig.pdb"]
CONTROLS_DIR = HERE / "controls"

# config name -> overrides of (dotted) keys in blastermaster_config.yaml
CONFIGS = {
    "default": {},
    "no_thin_spheres": {
        "dock_files_generation.thin_spheres_elec.use": False,
        "dock_files_generation.thin_spheres_desolv.use": False,
    },
}


def is_controlled(file_name):
    """Whether a file produced by a job is kept in the controls.

    Excluded: the untrimmed phi map (29 MB; trim.electrostatics.phi is derived from it)
    and the .dx visualization files (pure-Python conversions of the grids).
    """
    return file_name != "qnifft.electrostatics.phi" and not file_name.endswith(".dx")


@contextmanager
def chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def new_job(root, config_name):
    """Create a blastermaster job in `root`/job from the fixtures, with `config_name` applied."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    for f in FIXTURE_FILES:
        shutil.copy(f, root)
    with chdir(root):  # `new` copies input files from the current directory
        Blastermaster().new("job")

    job_dir = root / "job"
    config_path = job_dir / Blastermaster.CONFIG_FILE_NAME
    config = yaml.safe_load(config_path.read_text())
    for dotted_key, value in CONFIGS[config_name].items():
        *parents, key = dotted_key.split(".")
        d = config
        for p in parents:
            d = d[p]
        d[key] = value
    config_path.write_text(yaml.dump(config))
    return job_dir


def run_job(root, config_name):
    """Run a whole job; return {file name: path} of the controlled files it produced."""
    job_dir = new_job(root, config_name)
    working_dir = job_dir / Blastermaster.WORKING_DIR_NAME
    initial_files = set(os.listdir(working_dir))

    Blastermaster().run(str(job_dir))

    produced = {
        p.name: p
        for p in working_dir.iterdir()
        if p.is_file() and p.name not in initial_files and is_controlled(p.name)
    }
    produced[INDOCK_FILE_NAME] = job_dir / Blastermaster.DOCK_FILES_DIR_NAME / INDOCK_FILE_NAME
    return produced


def run_steps_in_isolation(root, config_name, control_dir):
    """Run each step of a job on control copies of its input files.

    Returns {file name: path} of the controlled files the steps produced, so that each step's
    output can be compared with the controls independently of the other steps.
    """
    job_dir = new_job(root, config_name)
    _, _, steps = load_steps(str(job_dir), str(job_dir / Blastermaster.CONFIG_FILE_NAME))
    outputs_dir = Path(root) / "outputs"
    outputs_dir.mkdir()
    produced = {}
    for step in steps:
        for infile in step.infiles:
            control = Path(control_dir) / infile.original_file_in_working_dir.name
            if control.exists():
                shutil.copy(control, infile.original_file_in_working_dir.path)
        step.run()
        for outfile in step.outfiles:
            name = outfile.original_file_in_working_dir.name
            if is_controlled(name):
                produced[name] = Path(shutil.copy(outfile.original_file_in_working_dir.path, outputs_dir / name))
    return produced
