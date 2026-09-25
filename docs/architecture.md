# Architecture

## Repository layout

```
pydock3/                  # the Python package
├── scripts.py            # the `pydock3` command
├── blastermaster/        # receptor + ligand → DOCK input files ("dockfiles")      → blastermaster.md
├── dockopt/              # optimize dockfile/INDOCK parameters by retrospective docking → other-tools.md
├── retrodock/            # retrospective docking of actives and decoys              → other-tools.md
├── criterion/            # enrichment metrics (ROC, normalized LogAUC, Bonferroni)
├── docking/              # DOCK 3.8 (git submodule dock3/, with the dock64 binary) and rundock.bash
├── lsd/                  # stub
├── jobs.py, job_schedulers.py   # array docking jobs on Slurm / SGE
├── config.py             # Parameter, ParametersConfiguration (YAML + yamale schema)
├── files.py              # File/Dir and file formats (INDOCK, OUTDOCK, mol2, db2 tarballs, ...)
├── util.py               # Script base class, logging, system_call, helpers
└── top_poses.py          # merge docked poses (not on the command line)
native/                   # sources of blastermaster's compiled programs  → native-programs.md
CMakeLists.txt            # builds them                                     → build-and-ci.md
pyproject.toml            # package metadata and build configuration
tests/                    # fixtures (rec.pdb, xtal-lig.pdb) and blastermaster's tests → testing.md
docs/                     # this documentation
.github/workflows/        # CI: wheels for each platform, tested
```

## The command line

`pydock3` (`scripts.py:main`) uses [python-fire](https://github.com/google/python-fire):
`pydock3 <script> [constructor args] - <command> [args]`. The script is one of `blastermaster`,
`dockopt`, `retrodock` (classes deriving from `util.Script`); the `-` separates the script's
(empty) constructor arguments from the command, and each command is a method:

```bash
pydock3 blastermaster - new              # Blastermaster.new(job_dir_path="blastermaster_job")
pydock3 blastermaster - run --help       # the arguments of Blastermaster.run
pydock3 dockopt - run slurm              # Dockopt.run(scheduler="slurm")
```

Every script follows the same pattern: `new` creates a job directory with a config file (from
the package's default config), copying input files from the current directory; `run` runs the
job from its directory.

## Core building blocks

- **Configuration** (`config.py`): YAML files validated with [yamale](https://github.com/23andMe/Yamale)
  schemas (`*_config_schema.yaml`), flattened into `Parameter(name, value)` objects with dotted
  names (`dock_files_generation.box_generation.margin`). Parameters hash deterministically, which
  dockopt uses to identify steps and results.
- **Files** (`files.py`): `File` and `Dir` wrap paths with existence checks, copying and
  resetting; `IndockFile` writes DOCK's `INDOCK`; readers for DOCK's outputs (`OutdockFile`),
  mol2, db2 and tarballs of ligands.
- **Blastermaster steps** (`blastermaster/util.py`): `BlasterStep` (a unit of work with declared
  input files, output files and parameters, run in its own directory and skipped when its outputs
  exist), `BlasterFiles` (all the files of a job), `program_path`. dockopt reuses them to build
  and run many variants of the pipeline ([blastermaster.md](blastermaster.md#how-a-step-runs)).
- **Jobs** (`jobs.py`, `job_schedulers.py`): array docking jobs and the Slurm/SGE interfaces
  ([other-tools.md](other-tools.md)).

## Dependencies

Python: numpy, scipy, pandas (data), networkx (dockopt's step graph), matplotlib, seaborn,
joypy, plotly (plots and reports), rdkit (SMILES checks), PyYAML/oyaml/yamale (configs), fire
(CLI), xmltodict (SGE's `qstat -xml`), timeout-decorator (jobs). python-dotenv, tornado and
pillow are declared in `pyproject.toml` but not imported by pydock3.

Compiled, bundled: the programs in `native/` (reduce, dms, filt, sphgen, chemgrid, qnifft,
solvmap), built at install time. Not bundled: DOCK itself (`dock64`, from the `dock3` submodule
or supplied by the user), only needed by dockopt and retrodock.

## Platform support

| | Linux x86_64 / aarch64 | macOS arm64 / x86_64 | Windows x86_64 |
|---|---|---|---|
| blastermaster | ✓ | ✓ | ✓ |
| dockopt, retrodock | ✓ (x86_64, with Slurm or SGE) | – | – |

How blastermaster got there: [cross-platform.md](cross-platform.md).
