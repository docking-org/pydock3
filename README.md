# pydock3

**pydock3** is a Python package that wraps the Fortran program UCSF DOCK (3.8 and later) and
provides tools to standardize and automate the steps of a molecular docking campaign: preparing a receptor for docking, and optimizing and evaluating docking set-ups. It is a
successor to *DOCK Blaster* (2009) and *blastermaster* (DOCK 3.7, 2012).

| Tool | What it does | Runs on |
|---|---|---|
| **blastermaster** | receptor + bound ligand → everything DOCK needs to dock (scoring grids, matching spheres, protonated receptor, `INDOCK`) | Linux (x86_64, aarch64), macOS (Apple silicon, Intel), Windows (x86_64) |
| **dockopt** | tries many blastermaster/DOCK parameter combinations and picks the best by docking known actives and decoys | Linux cluster (Slurm or SGE) |
| **retrodock** | docks actives and decoys with one set-up and reports enrichment | Linux cluster (Slurm or SGE) |

Documentation: [docs/](docs/README.md) (for maintainers and contributors), and the
[pydock3 wiki page](https://wiki.docking.org/index.php/Pydock3).

## Install

```bash
pip install "git+https://github.com/docking-org/pydock3.git@cross-platform-blastermaster"
```

This compiles blastermaster's programs, so it needs a C/C++ compiler and gfortran
([how to get them](docs/build-and-ci.md#installing-from-source)). Or install a prebuilt wheel,
which needs no compiler: download one from the latest run of the
[wheels workflow](https://github.com/docking-org/pydock3/actions/workflows/wheels.yml)
(*Artifacts*; needs a GitHub login), then `pip install pydock3-*.whl`. Wheels exist for Linux x86_64/aarch64
(glibc ≥ 2.28), macOS arm64 (≥ 14) and x86_64 (≥ 15), and Windows x86_64. Python 3.8–3.10.

## Quick start: prepare a receptor for docking

```bash
cd my_target/                    # contains rec.pdb (receptor) and xtal-lig.pdb (bound ligand)
pydock3 blastermaster - new      # creates blastermaster_job/
cd blastermaster_job
pydock3 blastermaster - run      # ~10 minutes for the test receptor
```

The files DOCK needs are then in `blastermaster_job/dockfiles/`, with `INDOCK` pointing at
them; `visualization/` has the grids and spheres for a molecular viewer. Parameters are in
`blastermaster_config.yaml`; to use your own version of any intermediate file (e.g. a protonated
receptor `rec.crg.pdb`), put it next to `rec.pdb` before `new`. See
[docs/blastermaster.md](docs/blastermaster.md).

dockopt and retrodock work the same way (`pydock3 dockopt - new`, `pydock3 dockopt - run slurm`),
with `actives.tgz` and `decoys.tgz`; see [docs/other-tools.md](docs/other-tools.md).

## What's new: blastermaster on any OS and CPU

Blastermaster used to run only on Linux x86-64: it drove prebuilt 32-bit Linux binaries through
shell commands, Perl and csh. Now:

- **No shell, Perl or csh.** Programs are run directly from Python; the Perl scripts and two small
  Fortran programs are ported to Python, byte-for-byte faithful.
- **The programs are built from source** (`native/`) for each platform, with a few portability
  fixes; reduce is vendored from its upstream source.
- **The same results everywhere.** The build avoids floating-point differences between CPUs
  and compilers, and dms uses correctly rounded math functions; every platform reproduces the
  same control files (within small tolerances for three of the grids).
- **Tested.** The tests compare every file blastermaster produces with control files, end to end
  and step by step; CI builds and tests a wheel for each platform.
- **Compared to the old binaries**, the results are the same up to floating-point rounding,
  except that sphgen makes a few borderline sphere choices differently (2 of 45 matching spheres
  change on the test receptor).

The whole story, and why each choice was made: [docs/cross-platform.md](docs/cross-platform.md).

## Repository

```
pydock3/          the package: blastermaster/, dockopt/, retrodock/, criterion/, docking/, ...
native/           sources of the compiled programs (reduce, dms, sphgen, chemgrid, qnifft, solvmap, filt)
tests/            test receptor and ligand, blastermaster's tests and control files
docs/             documentation
CMakeLists.txt, pyproject.toml, .github/workflows/wheels.yml   build and CI
```

## Development

```bash
pip install -e ".[dev]"          # editable install (rebuild after changing native/: run it again)
pytest tests/blastermaster       # ~30 minutes
```

Start with [docs/README.md](docs/README.md); common changes are in
[docs/maintenance.md](docs/maintenance.md).

## License

pydock3 is GPL-3.0-or-later ([LICENSE](LICENSE)). The bundled programs in `native/` keep their own
licenses ([docs/native-programs.md](docs/native-programs.md#licenses)).
