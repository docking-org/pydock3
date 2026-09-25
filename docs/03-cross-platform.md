# Making blastermaster OS- and CPU-agnostic

`pip install pydock3` now gives a working blastermaster on:

| Platform | Wheel tag | Toolchain |
|---|---|---|
| Linux x86_64, aarch64 | `manylinux_2_28` | gcc/gfortran 14 (manylinux image) |
| macOS arm64 | `macosx_14_0_arm64` | Homebrew gcc/gfortran 14 |
| macOS x86_64 | `macosx_15_0_x86_64` | Homebrew gcc/gfortran 14 |
| Windows x86_64 | `win_amd64` | MinGW-w64 (UCRT64) gcc/gfortran |

and every platform's wheel is tested in CI against the same control files. This page explains
what stood in the way and how each problem was solved, so that the reasoning survives.

## What made it Linux-only

1. **Shell glue.** Every step ran commands through `subprocess.run(..., shell=True)`: `sed -i`,
   `grep`, `cat`, redirections, pipes, `rm -rf`, Perl scripts run by their shebang, and a csh
   wrapper with a here-document.
2. **Prebuilt binaries.** The programs were shipped as Linux ELF binaries, mostly 32-bit static
   PGI builds, which cannot run on ARM, macOS or Windows. dms also found its helper through
   `/proc`, and forked it.
3. **Packaging.** poetry built a pure-Python (`py3-none-any`) wheel, so pip happily installed
   those Linux binaries anywhere.

## How it was done

The order mattered: **first make it measurable, then change it.**

1. **Controls.** Before any change, the original pipeline was run on a test receptor for two
   configurations and all its outputs saved (`tests/blastermaster/controls/legacy`; two runs were
   byte-identical, so the baseline is deterministic).
2. **Remove the scripting**, keeping the old binaries. With the same binaries, the output had to
   stay **byte-identical**, and it did (all 50 files), which isolates this step from any
   numerical change.
3. **Build the programs from source**, for every platform. Now differences are numerical only,
   and were reviewed program by program (below).
4. **Wheels and CI** for each platform, tested against controls made with the new build
   (`controls/rebuilt`).

### 1. No shell, no Perl, no csh

- `BlasterStep.run_program()` runs programs directly (`subprocess.run` with an argument list,
  no shell), with stdin/stdout redirected to files when a program needs it, and checks exit codes
  (before, failures only showed up as missing files later).
- The `sed`/`grep`/`cat` text processing became a few lines of Python in the steps, matching the
  tools' exact behaviour (e.g. `sed` works line by line, `grep` adds a final newline).
- The Perl scripts (makebox, makespheres1, makespheres3) and the tiny Fortran programs pdbtosph
  and showsphere became Python, ported line by line and checked byte for byte against the
  originals on two receptors over a sweep of parameters (170 cases). Details in
  [02-native-programs.md](02-native-programs.md#python-ports-of-former-scripts-and-small-programs).
- `rm -rf` became `shutil.rmtree`.

### 2. Paths and files

- Programs are found with `program_path(name)`, which looks for `bin/<name>` (`.exe` on Windows)
  in the package directories (there are two in an editable install).
- `INDOCK` refers to files with `/` (DOCK reads it; Windows accepts `/` too).
- Files blastermaster writes use LF line endings. Programs compiled for Windows write CRLF in
  text files; the text dockfiles are converted to LF when copied to `dockfiles/`, so they are the
  same on every platform (`vdw.vdw` and `.phi` are binary and never touched).
- Step directories are named `<step>_<first outfile>` instead of listing every outfile
  (~145 characters), because on Windows a process's working directory is limited to ~258
  characters.
- The Fortran programs truncate file names to 60–80 characters, so they are always given bare
  names and run in the step directory.

### 3. Building the programs

Sources live in `native/` and are built by CMake through scikit-build-core, with the same flags
everywhere (see [04-build-and-ci.md](04-build-and-ci.md)). The sources needed a few fixes to compile
and run as 64-bit programs on all platforms: a missing header (chemgrid), 32-bit pointer
assumptions and non-standard I/O (solvmap), and Linux/Unix-only process handling (dms). They are
listed in [02-native-programs.md](02-native-programs.md). Two flags matter for behaviour rather than
portability:

- `-fno-automatic`: the Fortran was written for compilers (PGI) that give local variables static
  storage, zero-initialized and kept between calls. solvmap relies on it (uninitialized local
  pointers passed to `realloc`) and crashed without it.
- `-fconvert=big-endian` for chemgrid and qnifft, whose grids DOCK reads as big-endian.

### 4. Wheels

Each wheel holds executables, not Python extensions, so it is tagged `py3-none-<platform>`: one
wheel per platform, for any Python 3. Compiler runtime libraries are linked statically (fully
static on Windows); macOS wheels bundle Homebrew's `libquadmath`. See
[04-build-and-ci.md](04-build-and-ci.md).

## Floating-point reproducibility

The same source can compute slightly different numbers on different machines. Blastermaster
makes discrete decisions from computed values (is a point on the surface? is a sphere too
close?), so a difference in the last bit can occasionally flip a decision, and flips propagate
downstream. The causes, and what the build does about each:

| Cause | Effect | What we do |
|---|---|---|
| **x87 80-bit arithmetic, aggressive optimization** (the old PGI `-tp p6 -fast` builds) | intermediate results rounded differently | nothing to do going forward: all current targets use IEEE single/double arithmetic (SSE, NEON) |
| **Fused multiply-add** (ARM CPUs, and compilers contracting `a*b+c` by default) | `a*b+c` rounded once instead of twice | `-ffp-contract=off` everywhere |
| **Vectorized math functions** (gfortran calls glibc's `libmvec` for vectorized loops, e.g. `expf` in qnifft) | less accurate than the scalar functions, x86 Linux only | `-fno-tree-vectorize` for the Fortran code |
| **Math library rounding** (`sin`, `exp`, `acos`, ... are not exactly rounded, and each OS's library differs in the last bit) | e.g. ~0.1% of dms's surface points typed differently on macOS/Windows | dmsd uses correctly rounded functions (CORE-MATH), identical everywhere; for qnifft and solvmap the differences are continuous and within the tests' tolerances |
| `-ffast-math` and friends | reordered arithmetic | never used |

`+ − × ÷ √` are exactly rounded by IEEE 754, so with these measures the only remaining source of
difference is the transcendental functions of qnifft (`exp`, `sin`, `cos`) and solvmap (`sin`,
`cos`). In practice: a local gfortran 8.5 build and the manylinux gfortran 14.2 build produce
bit-identical output, and in CI every platform's step-by-step outputs match the Linux controls
within the tolerances in [05-testing.md](05-testing.md).

## Differences from the legacy binaries

The rebuilt programs match the legacy ones step by step (each program given the legacy inputs):

| Program | Difference |
|---|---|
| reduce, filt, dms, and all Python steps | none (byte-identical) |
| chemgrid (`vdw.vdw`) | float32 rounding: relative ≤ 1e-4, median 8e-8 |
| qnifft (`trim.electrostatics.phi`) | ≤ 0.011 kT/e (values span ±1000) |
| solvmap (`ligand.desolv.*`) | 1 value in 210,000 differs, by 1 in the last printed digit |
| sphgen (`all_spheres.sph`) | coordinates within 1e-4 Å, but 0.5% of borderline spheres (43 of 7753) are kept or dropped differently |

End to end, sphgen's differences change 2 of the 45 matching spheres (one moves 0.5 Å, one is
replaced); with thin spheres off, they also change 27 of 120 low-dielectric spheres and so the
electrostatic potential near them (median 0.01 kT/e, ≤ 41 kT/e at a few grid points). No gfortran
setting (including x87 arithmetic) reproduces PGI's choices, and neither build is more correct
than the other. The source-built programs were adopted as the reference (`controls/rebuilt`);
`test_rebuilt_controls_match_legacy` keeps checking every file not downstream of sphgen against
the legacy controls.

## Platform notes

**Linux.** Built in the `manylinux_2_28` image (glibc ≥ 2.28: RHEL 8+, Ubuntu 18.10+). The
compiler runtimes (libgfortran, libstdc++, libgcc) are linked statically; auditwheel checks that
nothing else is needed.

**macOS.** gcc, g++ and gfortran 14 from Homebrew (installed by `fortran-lang/setup-fortran`,
which also sets `CC`/`CXX`/`FC`). Homebrew's runtime is built
for the runner's macOS version, so wheels require macOS 14 (arm64) / 15 (x86_64); the deployment
target must be set through `CIBW_ENVIRONMENT_MACOS` (cibuildwheel overrides a plain environment
variable). delocate copies `libquadmath` into the wheel. GitHub's Intel macOS runners are 2–3×
slower than the others for every step (a known runner problem); the tests take ~70 minutes there.

**Windows.** MinGW-w64 gcc/g++/gfortran (MSYS2 UCRT64), linked fully static so the programs need
only system DLLs (CI checks). Things that needed care:

- no `fork`: dms starts dmsd with `_spawnl`, handing it the pipes as stdin/stdout (with
  `_O_NOINHERIT` on our ends, or dmsd never sees end-of-file);
- stdio defaults to text mode, which would corrupt dms's binary protocol;
- C runtime differences (`creat` rejects Unix modes; no `caddr_t`, `sys/uio.h`, `roundeven`);
- CRLF line endings in the programs' text output (see above);
- cibuildwheel runs some commands through `cmd.exe`, where `<` in a requirement like
  `pytest<8` is a redirection.

## Not done (yet)

- **Windows on ARM**: no gfortran for it (LLVM flang would be needed), and `rdkit` has no
  win_arm64 wheels. Users can install the x64 wheel under Windows' x64 emulation.
- **musllinux** (Alpine): would be easy to add to CI.
- **dockopt, retrodock and DOCK itself** still need Linux, a job scheduler, and the `dock64`
  binary (see [06-other-tools.md](06-other-tools.md)).
- **Python ≥ 3.11**: `requires-python` and the dependency pins (e.g. `rdkit-pypi`, pandas 1.x) are
  unchanged from before.
