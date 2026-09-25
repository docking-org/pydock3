# dockopt, retrodock and the job schedulers

Besides blastermaster, pydock3 optimizes and evaluates docking set-ups by **retrospective
docking**: dock known actives and property-matched decoys, and measure how well the actives rank
(enrichment). These tools run DOCK itself, on a cluster, and are **still Linux-only**: they
haven't been part of the cross-platform port (see the end of this page).

| Module | What it does |
|---|---|
| `dockopt/` | searches over dockfile-generation and INDOCK parameters: builds each variant (with blastermaster's steps), docks actives and decoys with it, ranks the variants by enrichment |
| `retrodock/` | docks actives and decoys with one given set of dockfiles and reports enrichment |
| `criterion/` | enrichment metrics: ROC, normalized LogAUC, and Bonferroni-corrected significance thresholds (precomputed tables for 1–100 actives) |
| `jobs.py`, `job_schedulers.py` | array docking jobs on Slurm or SGE |
| `docking/` | DOCK 3.8 (`dock3/`, a git submodule with its Fortran sources and a prebuilt `dock64`) and `rundock.bash`, the script each array task runs |
| `lsd/` | a stub (config only; no command) |
| `top_poses.py` | merges docked poses into the top N by energy (not wired to the command line) |

## Requirements

- Linux with a **Slurm** or **SGE** cluster and a filesystem shared with the compute nodes.
- Environment variables:

  | Variable | Needed for |
  |---|---|
  | `TMPDIR` | scratch space on the compute nodes (rundock.bash stages files there) |
  | `SBATCH_EXEC`, `SQUEUE_EXEC` (absolute paths) | Slurm |
  | `QSUB_EXEC`, `QSTAT_EXEC` (absolute paths) | SGE |
  | `SLURM_SETTINGS`, `SGE_SETTINGS` (optional) | a script sourced before submitting |

- `dock64`: from the `dock3` submodule (`git submodule update --init`), an x86-64 Linux binary;
  or give your own with `custom_dock_executable` (dockopt) / `--custom_dock_executable`
  (retrodock). CI wheels don't include it (CI doesn't check out submodules).
- Inputs: `actives.tgz` and `decoys.tgz`, tarballs of DOCK `.db2`/`.db2.gz` ligand files.

## retrodock

```bash
cd my_target/        # dockfiles/, INDOCK, actives.tgz, decoys.tgz
pydock3 retrodock - new            # creates retrodock_job/
cd retrodock_job
pydock3 retrodock - run slurm      # or sge
```

It submits one array job for the actives and one for the decoys, waits, reads their `OUTDOCK`
files, sorts molecules by energy (ties rank decoys first) and writes the enrichment (normalized
LogAUC), a ROC plot, and energy and charge distributions.

## dockopt

```bash
cd my_target/        # rec.pdb, xtal-lig.pdb, actives.tgz, decoys.tgz
pydock3 dockopt - new              # creates dockopt_job/ with dockopt_config.yaml
cd dockopt_job                     # (dockopt must be run from inside its job directory)
pydock3 dockopt - run slurm        # or sge; see `pydock3 dockopt - run --help` for options
```

**Configuration** (`dockopt_config.yaml`): a `pipeline` of components, each a **step** (one round
of variants) or a **sequence** (steps repeated for `num_iterations`, stopping after
`max_iterations_with_no_improvement`). A step's `parameters` cover `dock_files_generation` (the
blastermaster parameters), `dock_files_modification` (random perturbation of the matching
spheres), `indock_file_generation` and `custom_dock_executable`. **Lists expand as a cartesian
product**, so every combination is tried. `^` means "the value from the previous component's best
result", and `{reference_value, arguments, operator}` builds a list relative to a value. Each
component keeps its `top_n` configurations by its `criterion` for the next.

**How it runs** (`dockopt.py`):

1. For each combination of dockfile parameters, `get_blaster_steps()` gives blastermaster's steps.
   All of them go into one graph (networkx) whose nodes are files, parameters and steps,
   identified by hashes, so steps shared by several variants run once.
2. The graph is run level by level. Steps marked `dockopt_submit_to_scheduler` (the desolvation
   grids, by default) are **pickled** and submitted as cluster jobs, which unpickle and run them;
   the others run locally. BlasterSteps must therefore stay picklable.
3. Each docking configuration (dockfiles × INDOCK parameters × DOCK executable) is docked against
   the actives and decoys as array jobs (`array_job_specs/*.txt`: one line per configuration);
   failed tasks are resubmitted.
4. The criterion is computed for each; `results.csv`, an HTML report (plotly, loaded from a CDN)
   and the best configurations (`best_retrodock_jobs/`, symlinks) are written.

**Job directory:** `dockopt_config.yaml`, inputs, `results.csv`, `report.html`,
`best_retrodock_jobs/`, `actives/` and `decoys/` (extracted), and one directory per component
(`1_step/`, or `2_seq/1_iter/1_step/` for sequences) with its `working/` (dockfiles with `_N`
suffixes, `INDOCK_N`, step directories) and `retrodock_jobs/`.

## Job schedulers (`jobs.py`, `job_schedulers.py`, `docking/rundock.bash`)

- `ArrayDockingJob` submits array jobs through `SlurmJobScheduler` (`sbatch --array`,
  `squeue`) or `SGEJobScheduler` (`qsub -t`, `qstat -xml`). A task is done when its `OUTDOCK.0`
  exists, failed when it has none and isn't queued.
- Each task runs `rundock.bash` (bash): finds its configuration line from the array task ID,
  stages the inputs in `$TMPDIR`, runs DOCK (forwarding SIGUSR1 so DOCK can write a restart
  file before a timeout), and copies back `OUTDOCK.N` and `test.mol2.gz.N`.
- Commands go through `pydock3.util.system_call` (`subprocess.run(shell=True, env=...)`); `env`
  replaces the environment, hence the absolute `*_EXEC` paths.
- `submit_single_step` writes a `submission.sh` that unpickles a BlasterStep and runs it with the
  same Python (`sys.executable`), which the compute nodes must be able to run.

## Linux-specific parts (for a future port)

- `system_call(shell=True)` everywhere in the schedulers; `source` (not POSIX `sh`); SGE options
  specific to UCSF's cluster (`-q !gpu.q`).
- `rundock.bash` (bash, awk, xargs, realpath, `kill -10` meaning SIGUSR1 on Linux).
- `jobs.py` uses `timeout_decorator` (SIGALRM: Unix only, main thread only).
- `dock64` is a Linux x86-64 PGI build; DOCK would need the same treatment as blastermaster's
  programs (or its Python reimplementation).
- `os.symlink` for the best results (needs privileges on Windows); `top_poses.py` uses `find`.

## Known quirks

- `pydock3 SDIFile ...` doesn't work: script names are matched lowercased in one place and not in
  another (`scripts.py`).
- With SGE, submitting pickled steps fails: the job name is an integer (`id(step)`), which the
  SGE code expects to be a string.
- dockopt and retrodock extract `actives/` and `decoys/` relative to the current directory: run
  them from inside the job directory.
- The config schema allows only `+` and `*` as operators, while `parameters.py` implements
  `+ - * /`.
- `rundock.bash`'s check for required variables never triggers (`[[ -z var ]]` tests the literal
  string).
- `util.CleanExit` logs an exception and swallows it, so a failed run still exits with status 0.
