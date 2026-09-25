# Blastermaster

Blastermaster turns a receptor structure (`rec.pdb`) and a bound ligand (`xtal-lig.pdb`) into
the files UCSF DOCK 3.8 docks against (the "dockfiles"): scoring grids, matching spheres, a
protonated receptor and an `INDOCK` parameter file. It is a pipeline of 24 steps (default configuration); 10 of them run bundled scientific
programs (see [native-programs.md](native-programs.md)), the rest are Python.

## Usage

```bash
cd my_target/                    # contains rec.pdb and xtal-lig.pdb
pydock3 blastermaster - new      # creates blastermaster_job/
cd blastermaster_job
# optional: edit blastermaster_config.yaml
pydock3 blastermaster - run      # ~10 min
```

(`-` separates the script name from its command; see [architecture.md](architecture.md#the-command-line).)

`new` creates the job directory:

```
blastermaster_job/
├── blastermaster_config.yaml   # parameters (from default_blastermaster_config.yaml)
├── working/                    # all files: inputs, defaults, intermediates, one dir per step
├── dockfiles/                  # the output: what DOCK needs
└── visualization/              # .dx grids and .pdb spheres to look at in a viewer
```

and copies into `working/`: `rec.pdb` and `xtal-lig.pdb` from the current directory, and the
parameter files from `pydock3/blastermaster/defaults/` (charges, radii, reduce's dictionary, ...).
**Any file in the current directory named like a blastermaster file is copied instead of the
default, and a step whose outputs already exist is skipped.** That is how you customize a run:
bring your own `rec.crg.pdb` (protonated receptor), `box`, `matching_spheres.sph`, parameter files...

`run` executes the steps whose outputs don't exist yet, then copies the dockfiles to `dockfiles/`
and writes `dockfiles/INDOCK`. Re-running is incremental: delete a file from `working/` and
the step that makes it runs again. Steps downstream are *not* re-run unless you delete their
outputs too.

## Outputs (`dockfiles/`)

| File | What it is | Made by |
|---|---|---|
| `rec.crg.pdb` | Receptor with hydrogens (polar only) and DOCK residue names | reduce + Python |
| `matching_spheres.sph` | Points ligand atoms are matched to, to place ligands in the site | makespheres (Python) |
| `vdw.vdw`, `vdw.bmp` | van der Waals scoring grid (repulsive/attractive terms) and bump map | chemgrid |
| `trim.electrostatics.phi`, `phi.size` | Electrostatic potential (Poisson–Boltzmann), trimmed to the box | qnifft + Python |
| `ligand.desolv.heavy`, `ligand.desolv.hydrogen` | Ligand desolvation grids (heavy atoms, hydrogens) | solvmap |
| `vdw.parms.amb.mindock` | van der Waals parameters (copied) | default file |
| `xtal-lig.pdb` | The ligand (copied) | input |
| `INDOCK` | DOCK's parameter file, pointing at the files above as `../dockfiles/<name>` | Python |

`vdw.vdw` and the `.phi` files are **big-endian** Fortran unformatted files (DOCK reads them so);
the others are text with LF line endings on every platform.

`visualization/` has the grids as OpenDX files (`vdw*.dx`, `ligdesolv.dx`,
`trim.electrostatics.dx`) and the matching spheres as `matching_spheres.pdb`.

## The pipeline

Steps, in order, for the default configuration. File names are the files in `working/`.

| # | Step (class, in `steps/`) | Runs | Inputs → outputs |
|---|---|---|---|
| 1 | `ReceptorMostOccupiedResiduesRenamingStep` | Python (`pdb.py`) | `rec.pdb` → `rec.most_occ_renamed.pdb` (keeps the most occupied alternate locations, clears alt/insertion codes, fixes chain IDs) |
| 2 | `ReceptorProtonationStep` | **reduce**, Python | → `rec.crg.pdb` (adds hydrogens, drops non-polar ones, renames HIS/CYS by protonation state) |
| 3 | `LigandHetatmRenamingStep` | Python | `xtal-lig.pdb` → `xtal-lig.hetatm_renamed.pdb` |
| 4 | `BindingSiteResiduesSelectionStep` | **filt** | receptor, ligand, `filt.params` → `rec.site` (residues near the ligand) |
| 5 | `MolecularSurfaceGenerationStep` | **dms** | `rec.crg.pdb`, `rec.site`, `radii` → `rec.ms` (molecular surface of the site) |
| 6 | `BindingSiteSpheresGenerationStep` | **sphgen** | `rec.ms` → `all_spheres.sph` (spheres filling the pocket) |
| 7 | `LigandPDBToSpheresConversionStep` | Python (`sphere_files.py`) | ligand → `xtal-lig.match.sph` (a sphere per heavy atom) |
| 8 | `MatchingSpheresGenerationStep` | Python (`makespheres.py`) | → `matching_spheres.sph` |
| 9 | thin spheres for electrostatics: `MolecularSurfaceGenerationStep`, `ThinSpheresGenerationStep`, `CloseSpheresGenerationStep`, `SpheresToPDBConversionStep` | **dms**, Python | → `rec.ts_elec.ms` → `thin_spheres_elec.sph` → `.sph.close` → `.sph.close.pdb` |
| 10 | thin spheres for desolvation (same four steps) | **dms**, Python | → `rec.ts_desolv.ms` → … → `thin_spheres_desolv.close.pdb` |
| 11 | `BoxGenerationStep` | Python (`makebox.py`) | ligand spheres → `box` (the grid box: ligand + 10 Å margin) |
| 12 | `ReceptorTransformationForElectrostatics` | Python | receptor + thin spheres → `receptor.crg.lowdielectric.pdb` |
| 13 | `ElectrostaticsGridGenerationStepYesThinSpheres` | **qnifft**, Python (`phi.py`) | → `qnifft.electrostatics.phi` → `trim.electrostatics.phi`, `phi.size`, `qnifft.atm` |
| 14 | `VDWScoringGridGenerationStep` | **chemgrid** | → `vdw.vdw`, `vdw.bmp` |
| 15 | `ReceptorTransformationForLigandDesolvationYesThinSpheres` | Python | receptor + thin spheres (as atom type X) → `rec.crg.lds.pdb` |
| 16 | `Hydrogen…`/`HeavyAtomLigandDesolvationScoringGridGenerationStep` | **solvmap** ×2 | → `ligand.desolv.hydrogen`, `ligand.desolv.heavy` |
| 17 | `VisualizationStep` | Python (`programs/visualization`) | grids → `.dx` files, `matching_spheres.pdb` |

Variants selected by the config:

- `thin_spheres_elec.use: false` replaces step 9 with `LowDielectricSpheresSelectionStep`
  (`makespheres.py`: pocket spheres from sphgen → `lowdielectric.sph`) and uses the
  `…NoThinSpheres` electrostatics step.
- `thin_spheres_desolv.use: false` skips step 10 and uses the receptor as is for desolvation.
- `covalent.use: true` adds `ChargedReceptorDeprotonationStep` after step 2 (removes the
  covalent residue's protons → `rec.crg.deprotonated.pdb`) and makes the matching spheres from
  that residue's atoms instead.

**Thin spheres.** A thin layer of spheres on the receptor's molecular surface, kept where they
are close to the ligand. For electrostatics they extend the receptor's low-dielectric region
into the site; for desolvation they are added as extra atoms. Their size and placement come from
the `thin_spheres_*` parameters.

The pipeline is assembled in `get_blaster_steps()` (`blastermaster.py`); `load_steps(job_dir,
config)` returns the steps of a job, which the tests use too.

## Configuration

`blastermaster_config.yaml` has two sections, validated against
`blastermaster_config_schema.yaml` (yamale):

- `dock_files_generation`: how the dockfiles are made. The file written by `new` has the main
  options (`receptor_protonation.reduce_options`, `thin_spheres_elec/desolv.*`, `covalent.*`).
  Optional sections override defaults defined in the step classes:

  | Key | Default | Used by |
  |---|---|---|
  | `molecular_surface_density` | 5.0 points/Å² | dms (main surface) |
  | `box_generation.margin` | 10.0 Å | makebox |
  | `matching_spheres_generation.{max_num_spheres, gridsize, tooclose}` | 45, 1.5, 0.8 | makespheres |
  | `low_dielectric_sphere_selection.min_num_spheres` | 25 | makespheres (no thin spheres) |
  | `electrostatics_grid_gen.{grid_size, use_receptor_box}` | 193, false | qnifft |
  | `vdw_grid_gen.grid_spacing` | 0.2 Å | chemgrid |
  | `desolv_grid_gen.{other_radius, probe_radius, hydrogen_radius, heavy_radius}` | 1.0, 1.4, 1.0, 1.8 Å | solvmap |

- `indock_file_generation`: values written into `INDOCK` (DOCK's search and scoring settings).

## How a step runs

All steps derive from `BlasterStep` (`pydock3/blastermaster/util.py`):

- A step declares **infiles, outfiles and parameters**. Each file is a `BlasterFile` in
  `working/`; the step gets its own copy in a **step directory**
  `working/<step_name>_<first outfile>/`, which is where it runs.
- `run()` is wrapped by `BlasterStep.handle_run_func`: if all outfiles exist in `working/` the
  step is skipped; otherwise the step dir is reset, infiles are copied in, the step runs, and its
  outfiles are copied back to `working/`.
- Programs are run with `run_program([program_path(name), args...], stdin_file_path=...,
  stdout_file_path=...)`: no shell, the step dir as working directory, the command, output and
  exit code appended to the step's `log`, and an error if the exit code isn't expected.
  `program_path()` finds the executables in `pydock3/blastermaster/bin/`.
- Several programs read inputs from fixed file names (e.g. `INCHEM`, `INSEV`, `INSPH`) or ask
  questions on stdin (filt); steps write those files and rename infiles as the programs expect.
- The Fortran programs truncate file names (to 60–80 characters), which is why they are always
  given bare names, never paths.

To debug a step, look in its step dir: its inputs, outputs, the programs' own output files
(`OUTCHEM`, `OUTSEV`, ...) and `log`.

dockopt builds many variants of these steps (different parameters) and runs them through job
schedulers; step instances must therefore stay picklable (see [other-tools.md](other-tools.md)).

## Code map

```
pydock3/blastermaster/
├── blastermaster.py        # Blastermaster script (new/run), get_blaster_steps, load_steps
├── util.py                 # BlasterStep, BlasterFile(s), program_path, file-name tables
├── config.py               # BlastermasterParametersConfiguration (yaml + schema)
├── pdb.py, phi.py          # PDB processing; DelPhi phi map reading/trimming/writing
├── steps/                  # one module per step
├── programs/
│   ├── makebox.py, makespheres.py, perl_semantics.py   # Python ports of the Perl scripts
│   ├── sphere_files.py     # Python ports of pdbtosph and showsphere
│   ├── thinspheres/        # sphere/PDB helpers for the thin spheres steps
│   └── visualization/      # grids → OpenDX
├── defaults/               # default parameter files copied into jobs
├── bin/                    # the native programs; only in installed packages (built by CMake)
├── default_blastermaster_config.yaml
└── blastermaster_config_schema.yaml
```
