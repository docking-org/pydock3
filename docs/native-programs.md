# The native programs

Blastermaster runs seven compiled programs (plus dms's helper, dmsd). Their sources are in
`native/`, they are built by CMake (see [build-and-ci.md](build-and-ci.md)) and installed into
`pydock3/blastermaster/bin/`. This page says what each one does, where it comes from, and every
change made to it. Changes in the sources are marked with `pydock3` comments.

Until the cross-platform port they were shipped as prebuilt Linux binaries, mostly 32-bit builds
made with the PGI compilers; those binaries and the sources as they were are in the git history
(e.g. commit `a7bba56`, `pydock3/blastermaster/programs/`).

| Program | Language | Origin | Role in blastermaster |
|---|---|---|---|
| reduce | C++ | Richardson Lab, Duke (J. M. Word) | add hydrogens to the receptor |
| filt | Fortran 77 | UCSF DOCK | pick the binding-site residues |
| dms (+ dmsd) | C | UCSF Computer Graphics Lab | molecular surface |
| sphgen | Fortran 77 | UCSF DOCK (Kuntz lab) | spheres filling the pocket |
| chemgrid | Fortran 77 | UCSF DOCK (E. C. Meng) | van der Waals grid and bump map |
| qnifft 2.2 | Fortran 77 | Kim Sharp, U. Pennsylvania | electrostatic potential (Poisson–Boltzmann) |
| solvmap | Fortran 77 + C | Shoichet lab (UCSF) | ligand desolvation grids |

Fortran programs read their input from fixed file names or stdin and write output in fixed
formats; the steps that run them are listed in [blastermaster.md](blastermaster.md#the-pipeline).

## reduce

Adds hydrogens to the receptor and optimizes flips (Asn, Gln, His) and rotatable hydrogens;
run as `reduce -db reduce_wwPDB_het_dict.txt -HIS -FLIPs rec.most_occ_renamed.pdb` (the options
come from `receptor_protonation.reduce_options`).

- Source: `native/reduce/`, [rlabduke/reduce](https://github.com/rlabduke/reduce) at commit
  `b3aac6e`, which reports version `reduce.4.9.210817` like the binary it replaces and gives
  identical output. Only `libpdb/`, `toolclasses/`, `reduce_src/`, the license and README are
  vendored (not the Python bindings, the 66 MB het dictionary or build files; pydock3 has its own
  dictionary in `defaults/`).
- Unchanged. Built with our own `add_executable` (upstream's `reduce_src/CMakeLists.txt` needs
  Python and Boost); its `libpdb` and `toolclasses` CMake files are used as they are.
- Exit code 1 means some flip optimizations were abandoned (too many combinations); the output
  is still complete, so the step accepts it.
- Newer versions (4.15, now Apache-2.0) may protonate differently: upgrading means regenerating
  the controls ([maintenance.md](maintenance.md)).

## filt

Selects the receptor residues near the ligand (writes `rec.site`), answering its questions from
`filt.params` on stdin. `native/filt/filter1.f`, unchanged; needs `-ffixed-line-length-none`
(tab-formatted lines run past column 72).

## dms and dmsd

dms computes the molecular (Connolly/Richards) surface of the binding site: points with normals,
typed contact (`SC`), saddle (`SS`) or reentrant. It was written as "distributed MS": the client
`dms` sends work to compute servers `dmsd`, over TCP on other machines or, normally, to one
`dmsd` it starts itself and talks to through pipes (a binary protocol on the server's
stdin/stdout).

Changes (`native/dms/`), all for portability:

- The TCP mode is removed (`compute.c`); dms always starts one local dmsd. `get_server()` waits
  for it with a blocking read instead of `select()`, which doesn't work on Windows pipes.
- dmsd is found next to dms from `argv[0]` (`dms_paths.c`), instead of `/proc/<pid>/exe`,
  which only Linux has. pydock3 always runs dms by its full path.
- On Windows (`#ifdef _WIN32`), dmsd is started with `_spawnl` with the pipes on its
  stdin/stdout (there is no `fork`), and both use binary-mode stdio. `iovec.h` supplies
  `struct iovec` and `caddr_t`, which MinGW lacks.
- dmsd's daemon housekeeping is removed (`chdir /tmp`, a lock file, `nice`, `sbrk`).
- dms opens its log (`-g`) with `freopen` (Windows' `creat` rejects mode 0666).
- **dmsd uses correctly rounded `acos`, `atan2`, `sin`, `cos`** (`cr_math.h`, from
  `native/core-math/`). The platforms' math libraries round these differently in the last bit,
  which flipped ~0.1% of the surface points (their type, or a few points in or out) and so
  changed everything downstream. Correctly rounded results are the same everywhere, and here
  identical to glibc's, so the Linux output is unchanged. See
  [cross-platform.md](cross-platform.md#floating-point-reproducibility).
- The code is 1990s K&R C, compiled with `-std=gnu89` (implicit declarations are valid there;
  newer standards reject them).

## sphgen

Generates spheres filling the receptor's pockets from the molecular surface, and clusters them
(reads `INSPH`, writes `all_spheres.sph`). `native/sphgen/sphgen.f`, unchanged.

## chemgrid

Computes the van der Waals grid (the repulsive and attractive terms, at each grid point, of the
AMBER-like receptor parameters) and the bump map (reads `INCHEM`, writes `vdw.vdw`, `vdw.bmp`).

- `parmrec.h` was missing from the sources; it is taken from DOCK 3.6
  (`src/docktools/parmrec.h`, where `parmrec.f` is identical). Its common blocks have exactly the
  sizes in the old binary's symbol table.
- `aesthetics.f` (not called by chemgrid) is left out.
- Flags: `-fconvert=big-endian` (DOCK reads `vdw.vdw` as big-endian) and `-fwrapv` (the
  atom-type hash relies on integer overflow wrapping around).

## qnifft

Solves the (non-linear) Poisson–Boltzmann equation by finite differences for the receptor's
electrostatic potential (Kim Sharp's DelPhi descendant, version 2.2). Run as
`qnifft qnifft.parm`; writes a DelPhi phi map, which `phi.py` trims to the box.

- `native/qnifft/`: the 27 files of the v2.2 build, with `qdiffpar.h` for a 193³ grid (the grid
  size is compiled in and is part of the phi file format) and `lstmod.h` holding the original
  build stamp. Unchanged. Other versions, examples and utilities were not kept.
- Flags: `-ffixed-line-length-132`, `-fconvert=big-endian`.
- Its arrays are static (~720 MB of zero-initialized memory, allocated lazily by the OS).

## solvmap

Computes the ligand desolvation grids: at each grid point, how much a ligand atom (heavy or
hydrogen, by Born radius) would be desolvated by the receptor (reads `INSEV`, writes
`ligand.desolv.*`). The slowest program (~3 minutes per grid).

Changes (`native/solvmap/`):

- `memalloc.c` returns memory for Fortran (Cray) pointers; its result is declared `integer*8`
  (it was a 4-byte `integer`, which only works in 32-bit builds), its Fortran name mangling is
  fixed, and a two-argument `perror` call is corrected.
- `form='binary'` (a PGI extension) becomes `access='stream', form='unformatted'` for its
  (unused) `.plt` output; a format string missing its parentheses is fixed.
- `cube.f`: an unused allocation through a 4-byte "pointer" is removed (on 64-bit it passed
  garbage to `realloc`).
- `solvmap.f`: a write past the end of a `(2,3)` array (`edgesolv(3,i)`, in a branch taken only
  without a box file) is removed.
- Flag: `-fcray-pointer`.

## Python ports of former scripts and small programs

These used to be Perl scripts, a csh wrapper and two tiny Fortran programs. They are now Python,
in `pydock3/blastermaster/programs/`, ported line by line so their output is identical,
including the originals' quirks (commented where they matter):

| Module | Replaces | Does |
|---|---|---|
| `makebox.py` | `makebox.smallokay.pl` | the grid box around the ligand spheres |
| `makespheres.py` | `makespheres3.cli.pl` (`make_matching_spheres`), `makespheres1.cli.pl` (`make_low_dielectric_spheres`) | select spheres: near the ligand and receptor, H-bond geometry, thinned, continuous |
| `sphere_files.py` | Fortran `pdbtosph`; Fortran `showsphere` + `doshowsph.csh` | PDB → sphere file; sphere cluster → PDB |
| `perl_semantics.py` | — | how Perl turns whitespace-split tokens into numbers/strings |

`sphere_files.py` reads fixed columns into float32, like the Fortran REALs it replaces, so
values print the same (e.g. 134.296 prints as 134.29601 at 5 decimals). The ports were checked
byte for byte against the originals on two receptors over a sweep of their parameters (170
cases). One deliberate difference: makespheres3 visited polar spheres in Perl's random hash
order; the port uses sorted order (as makespheres1 did), which gives the same output here.

## CORE-MATH

`native/core-math/`: correctly rounded `acos`, `atan2`, `sin`, `cos` from the
[CORE-MATH project](https://core-math.gitlabpages.inria.fr/) (MIT license), used by dmsd. One
change: `acos.c` always uses its portable `roundeven` fallback (`__builtin_roundeven` can become
a call to `roundeven()`, which MinGW's C library lacks). See `native/core-math/README.md`.

## Licenses

Each directory keeps its original notices: reduce (`native/reduce/LICENSE.txt`, and
`libpdb/LICENSE.txt` for UC's libpdb), dms (UC license in its sources), CORE-MATH (MIT). The DOCK
Fortran programs and qnifft carry their authors' notices in their headers; qnifft's `README`
describes its terms (free for academic research).
