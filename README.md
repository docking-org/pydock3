# pydock3

**pydock3** is a Python package wrapping the Fortran program *UCSF DOCK* that provides tools to help standardize and automate the computational methods employed in molecular docking. It is a natural successor to *DOCK Blaster*, originally published in 2009, and *blastermaster*, part of the DOCK 3.7 release in 2012. 

Documentation: https://wiki.docking.org/index.php/Pydock3

## Installation

```bash
pip install .        # from a clone of this repository
pip install -e .     # for development
```

This compiles the programs blastermaster runs (in `native/`), so it needs a C/C++ compiler,
`gfortran`, and CMake (pip installs CMake and Ninja if they are missing). The wheels built by
CI (`.github/workflows/wheels.yml`) contain the compiled programs for Linux (x86_64, aarch64),
macOS (arm64, x86_64), and Windows (x86_64), and need no compilers.

## Tests

```bash
pytest tests/blastermaster  # ~30 min
```

These run blastermaster on `tests/rec.pdb` and `tests/xtal-lig.pdb` and compare its output with
the control files in `tests/blastermaster/controls` (see `tests/blastermaster/make_controls.py`).
