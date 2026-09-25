# Maintaining pydock3

How to make common changes safely, and what is known to be unfinished.

## Principles

- **Measure first.** Blastermaster's output is compared to control files in every test run, on
  every platform. Before changing behaviour, know what the tests say; after, explain every
  difference ([05-testing.md](05-testing.md)).
- **Keep results reproducible across platforms**: no `-ffast-math`, no FMA contraction, no
  vectorized math functions, and correctly rounded math where a program makes discrete decisions
  from transcendental functions ([03-cross-platform.md](03-cross-platform.md#floating-point-reproducibility)).
- **Ports are faithful.** The Python ports of the Perl scripts and the patched programs reproduce
  the originals, quirks included. Fix a quirk as a deliberate, separate change that updates the
  controls.
- **No shell.** Run programs with `BlasterStep.run_program` (argument lists, no `shell=True`);
  do text processing in Python.
- Mark changes to vendored or legacy code with a `pydock3` comment.

## Common changes

### Changing a blastermaster step

Steps are in `pydock3/blastermaster/steps/` and wired together in `get_blaster_steps()`
(`blastermaster.py`). A step declares its infiles, outfiles and parameters in `__init__` and does
its work in `run()` (decorated with `@BlasterStep.handle_run_func`), in its own step directory.
To add a file, add it to `BLASTER_FILE_IDENTIFIER_TO_PROPER_BLASTER_FILE_NAME_DICT`
(`util.py`); to add a parameter, give it a default in the step (`Parameter("dock_files_generation.…", default)`)
and add it to `blastermaster_config_schema.yaml`. See [01-blastermaster.md](01-blastermaster.md#how-a-step-runs).

### Changing a program

1. Edit the source in `native/<program>/`, with a `pydock3` comment.
2. Rebuild: `pip install -e .` (or `cmake --build` a build directory, see
   [04-build-and-ci.md](04-build-and-ci.md)).
3. Run the tests; update the controls if the output is meant to change.

### Adding a program

1. Put its sources in `native/<program>/` (with its license).
2. Add an `add_executable` to `CMakeLists.txt`, with the flags it needs, and add it to the
   `install(TARGETS ...)` list.
3. In the step, run it with `self.run_program([program_path("<program>"), ...])`.
4. Check it builds on Windows too: it must not need POSIX-only APIs (`fork`, `select` on pipes,
   `/proc`, ...). You can cross-compile with MinGW-w64 in a Debian container
   (`gcc-mingw-w64-x86-64`, `g++-mingw-w64-x86-64`, `gfortran-mingw-w64-x86-64`, and a CMake
   toolchain file setting `CMAKE_SYSTEM_NAME=Windows` and those compilers) and run the `.exe`
   files under Wine, before waiting for CI.

### Updating reduce

1. Copy `libpdb/`, `toolclasses/`, `reduce_src/`, `LICENSE*` and `README.md` of the new version
   into `native/reduce/`, leaving out `reduce_src/CMakeLists.txt`, `reduce_bpl.cpp`, `reduce.py`,
   `SConscript`s and `Makefile`s.
2. Update the source list of `reduce` in `CMakeLists.txt` if files were added or removed, and
   the commit in [02-native-programs.md](02-native-programs.md#reduce).
3. Expect different hydrogens: regenerate the controls and review the differences. Since 4.15
   reduce is Apache-2.0 licensed (update the license notes).

### Updating CORE-MATH

Copy `src/binary64/{acos,atan2,sin,cos}/*.c` and `atan2/tint.h` from
[CORE-MATH](https://gitlab.inria.fr/core-math/core-math), re-apply the `roundeven` change in
`acos.c` (see `native/core-math/README.md`), and update the commit there. Correctly rounded
results don't change, so neither should the controls.

### Changing dependencies or supported Python versions

Edit `[project]` in `pyproject.toml`. CI builds with Python 3.10 but the wheels work with any
Python 3 allowed by `requires-python`; installing them in CI's test environment is what checks
that the dependencies install on every platform.

### CI maintenance

- `macos-15-intel` is GitHub's last Intel macOS image (available until August 2027). When it
  goes, drop the Intel macOS wheel or cross-compile it on an arm64 runner.
- Keep actions (`actions/checkout`, `upload-artifact`, `setup-fortran`, `setup-msys2`,
  `cibuildwheel`) reasonably current; GitHub warns about deprecated ones in the run summary.
- To add musllinux wheels (Alpine), remove `skip = "*-musllinux_*"` from `pyproject.toml`
  and install gfortran in the musllinux image (e.g. `before-all = "apk add gfortran"` in a
  `[[tool.cibuildwheel.overrides]]` for `*-musllinux_*`).

## Releasing

Not set up yet. To publish to PyPI: bump `version` in `pyproject.toml`, tag the release, and add
a job to `wheels.yml` that downloads the wheel and sdist artifacts and uploads them with
`pypa/gh-action-pypi-publish` (PyPI "trusted publishing", no token needed) on tags.

## Known issues and unfinished work

- **dockopt, retrodock, DOCK**: still Linux-only (job schedulers through bash scripts, `dock64`
  binary, `timeout_decorator`); see [06-other-tools.md](06-other-tools.md).
- **Python ≥ 3.11 and dependency pins** (`rdkit-pypi`, pandas 1.x, numpy 1.x, ...) are unchanged
  from before the port.
- **Windows on ARM**, **musllinux**: no wheels (see [03-cross-platform.md](03-cross-platform.md#not-done-yet)).
- **Covalent docking** (`covalent.use: true`) is not covered by the tests.
- `programs/visualization/create_VDW_DX.py` reads `vdw.vdw` without skipping its Fortran record
  markers and assumes a 0.2 Å grid spacing, so the vdW `.dx` files are slightly off. Only the
  visualization is affected, not the dockfiles.
- `Blastermaster.new`/`run` create `visualization/` from the given path rather than the job
  directory's absolute path (harmless in practice).
- Inherited from the original makebox: a margin above ~20 Å never converges (even a single point
  gives a box with too many grid points).
