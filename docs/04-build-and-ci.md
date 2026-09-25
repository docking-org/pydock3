# Building, packaging and CI

pydock3 is a Python package with compiled programs inside. The Python code is packaged as is;
the programs in `native/` are compiled by CMake and installed into
`pydock3/blastermaster/bin/`. [scikit-build-core](https://scikit-build-core.readthedocs.io)
drives both from `pyproject.toml`.

```
pyproject.toml      # metadata, dependencies, scikit-build-core and cibuildwheel settings
CMakeLists.txt      # how to compile the programs in native/
native/             # their sources (see native-programs.md)
.github/workflows/wheels.yml   # CI: wheels for every platform, tested
```

## Installing from source

```bash
pip install .       # build and install
pip install -e .    # editable install, for development
```

This compiles the programs, so it needs:

| Platform | Compilers |
|---|---|
| Linux | gcc, g++, gfortran (e.g. `apt install build-essential gfortran`, `dnf install gcc-c++ gcc-gfortran`) |
| macOS | `brew install gcc` (provides `gfortran`; C/C++ can come from Xcode's command line tools) |
| Windows | MSYS2 UCRT64: `pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran`, then `CMAKE_GENERATOR=Ninja CC=gcc CXX=g++ FC=gfortran` with its `bin` on `PATH` |

CMake (≥ 3.15) and Ninja are installed by pip if missing.

**Editable installs.** `pip install -e .` builds the programs and installs them into
site-packages; the Python code is used from the repository. `program_path()` searches all the
package's directories, so it finds them. Changing Python code needs nothing; changing code in
`native/` or `CMakeLists.txt` needs `pip install -e .` again (the CMake build directory,
`build/`, is reused, so that's quick).

**Building without pip**, e.g. to work on a program:

```bash
cmake -S . -B build/dev -DCMAKE_BUILD_TYPE=Release
cmake --build build/dev -j 8      # executables in build/dev/
```

## CMakeLists.txt

One file builds everything:

- one `add_executable` per program (`filt`, `sphgen`, `chemgrid`, `solvmap`, `qnifft`, `dms`,
  `dmsd`, `reduce`), plus `core_math` (a static library linked into dmsd) and reduce's `pdb++`
  and `toolclasses` libraries (through their own CMake files);
- `install(TARGETS ... RUNTIME DESTINATION pydock3/blastermaster/bin)`.

Compiler flags, and why:

| Flag | Where | Why |
|---|---|---|
| `-ffp-contract=off` | everything | no fused multiply-adds, which ARM CPUs (and some compilers by default) use: same results on every CPU |
| `-fno-tree-vectorize` | Fortran | vectorized loops call glibc's `libmvec` math functions, less accurate than the scalar ones and x86-Linux-only |
| `-std=legacy` | Fortran | FORTRAN 77 with extensions |
| `-fno-automatic` | Fortran | static, zero-initialized local variables, as the code (written for PGI compilers) expects |
| `-fconvert=big-endian` | chemgrid, qnifft | DOCK reads their grids as big-endian |
| `-fwrapv` | chemgrid | its atom-type hash relies on integer overflow wrapping |
| `-ffixed-line-length-none`, `-132` | filt, qnifft | tab-formatted lines past column 72 |
| `-fcray-pointer` | solvmap | Cray pointers |
| `-std=gnu89 -w` | dms, dmsd | K&R C (implicit declarations), and its many warnings |
| `-static-libgcc -static-libgfortran -static-libstdc++` (`-static` on Windows) | when `PYDOCK3_STATIC_RUNTIME=ON` (CI) | wheels don't depend on the compiler's runtime libraries |

`PYDOCK3_STATIC_RUNTIME` is off by default because a local compiler may not have static
runtime libraries (e.g. the system gfortran of RHEL/Rocky 8).

The build type is Release (`-O3` for gcc). Don't add `-ffast-math` or similar: see
[03-cross-platform.md](03-cross-platform.md#floating-point-reproducibility).

## pyproject.toml

- `[project]`: standard (PEP 621) metadata, dependencies and the `pydock3` command. The
  dependency bounds are the former poetry `^` constraints, translated; `dev` extras hold the
  development tools (pytest, ...).
- `[build-system]`: `scikit-build-core`.
- `[tool.scikit-build]`:
  - `wheel.py-api = "py3"`: the wheel contains executables, not Python extensions, so it is
    tagged `py3-none-<platform>`: one wheel per platform works with any Python 3.
  - `build-dir = "build/{wheel_tag}"`: a persistent build directory (fast rebuilds).
  - `sdist.exclude`: the test controls (~220 MB) are left out of the sdist.
- `[tool.cibuildwheel]`: how CI builds and tests wheels (next section).

## CI: `.github/workflows/wheels.yml`

On every push (any branch) and pull request, except those changing only documentation, or on a
manual trigger, GitHub Actions builds and tests a wheel on each platform with [cibuildwheel](https://cibuildwheel.pypa.io), and builds the sdist.

| Runner | Wheel | Toolchain setup |
|---|---|---|
| `ubuntu-latest` | Linux x86_64 | none: cibuildwheel builds in the `manylinux_2_28` container (gcc-toolset 14) |
| `ubuntu-24.04-arm` | Linux aarch64 | same |
| `macos-14` | macOS arm64 | `fortran-lang/setup-fortran` (Homebrew gcc 14; sets `CC`/`CXX`/`FC`) |
| `macos-15-intel` | macOS x86_64 | same |
| `windows-latest` | Windows x86_64 | `msys2/setup-msys2` (UCRT64 gcc, gfortran), put on `PATH`; `CMAKE_GENERATOR=Ninja CC=gcc CXX=g++ FC=gfortran` |

For each, cibuildwheel:

1. builds the wheel once (`build = "cp310-*"`; `archs = ["native"]`, `skip = "*-musllinux_*"`),
   with `PYDOCK3_STATIC_RUNTIME=ON`;
2. repairs it: auditwheel on Linux (checks the programs need only glibc, retags
   `manylinux_2_28`), delocate on macOS (bundles Homebrew's `libquadmath`, retags for the
   deployment target); nothing on Windows, where a later step checks with `objdump` that the
   programs import no MinGW DLLs;
3. installs it with its dependencies in a clean environment and runs
   `pytest tests/blastermaster` against it (with live logging of each step and durations).

The macOS deployment target (14.0 on arm64, 15.0 on Intel, the runners' versions, which
Homebrew's runtime requires) is passed as `CIBW_ENVIRONMENT_MACOS`; cibuildwheel overrides a
plain `MACOSX_DEPLOYMENT_TARGET`.

Wheels and sdist are uploaded as workflow artifacts. Publishing to PyPI is not set up yet (see
[07-maintenance.md](07-maintenance.md#releasing)).

**Timing.** Building takes ~1 minute; the tests ~20–30 minutes, and ~70 minutes on the Intel
macOS runner, which is slow for everything.

**Reading CI results** with the GitHub CLI:

```bash
gh run list -b <branch>                        # runs
gh run view <run-id>                           # jobs
gh run view --job <job-id> --log-failed        # failed steps' output (once the run is done)
gh api repos/<owner>/<repo>/actions/jobs/<job-id>/logs   # a job's full log
```

Test failures also appear as annotations on the run (the `pytest-github-actions-annotate-failures`
plugin), with the list of differing files.

**Building a manylinux wheel locally** (needs Docker), exactly as CI does:

```bash
pipx run cibuildwheel --platform linux --output-dir wheelhouse .   # add CIBW_TEST_SKIP="*" to skip tests
```
