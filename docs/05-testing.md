# Testing

The tests (`tests/blastermaster/`) run blastermaster on a real receptor and compare everything
it produces with **control files**: the output of a known-good build. They are slow (~30 minutes)
but they cover every program and step, on every platform (CI runs them against each wheel).

```bash
pytest tests/blastermaster                            # all tests
pytest tests/blastermaster -k legacy                  # just the static legacy comparison (seconds)
pytest tests/blastermaster -o log_cli=true --log-cli-level=INFO   # show each step as it starts
```

## Fixtures and controls

- `tests/rec.pdb`, `tests/xtal-lig.pdb`: the neuropeptide S receptor (NPSR1) with a bound peptide
  (source in `tests/README.txt`).
- `tests/blastermaster/controls/<tier>/<config>/`: every file a job produces, plus `INDOCK` and
  the config used. Two configurations (`CONFIGS` in `harness.py`):
  - `default`: the default config (thin spheres on);
  - `no_thin_spheres`: thin spheres off, which exercises the other steps (low-dielectric spheres
    from makespheres1, the no-thin-spheres electrostatics and desolvation).

  Two tiers:
  - `rebuilt`: made with the programs built from source (the manylinux x86_64 wheel). **What
    the tests compare against.**
  - `legacy`: made with the original code and prebuilt programs, before the cross-platform
    port. Kept as the historical reference.

  Not kept (large, and derived from kept files): the untrimmed `qnifft.electrostatics.phi` and the
  `.dx` visualization files.

## The tests

`test_blastermaster.py`:

| Test | Runs | Checks |
|---|---|---|
| `test_end_to_end[config]` | a whole job (`new` + `run`) | all its files against `rebuilt` |
| `test_steps_in_isolation[config]` | each step once, with the **control** files as its inputs | each step's outputs against `rebuilt` |
| `test_rebuilt_controls_match_legacy[config]` | nothing (compares control files) | `rebuilt` against `legacy`, except files downstream of sphgen |

Step isolation is what makes failures readable: when a program behaves differently, only its own
outputs differ, instead of everything downstream. The end-to-end test checks that the steps fit
together (file names, order, the `INDOCK` file).

`harness.py` has the helpers: `new_job` (a job from the fixtures, with a config applied),
`run_job`, `run_steps_in_isolation` (uses `load_steps()` from `blastermaster.py`) and
`is_controlled` (which files are kept).

## How files are compared

`compare.py` compares files numerically rather than byte by byte:

- **Text files**: line by line, the numbers in each line must be within tolerance, and the rest of
  the line (with numbers and whitespace taken out) must be identical. So a number's width or the
  line endings may change, but not the structure or any text.
- **Binary grids** (`vdw.vdw`, `trim.electrostatics.phi`): big-endian Fortran records; the
  records holding float32 values are compared as numbers, the others (labels) exactly.

A failure lists each differing file, e.g.
`all_spheres.sph: 64202/64202 values differ (max abs diff 0); part 219: expected b'cluster#…', got b'########'`
(values shifted, because a line appeared or vanished).

## Tolerances

Everything must match exactly, except (`TOLERANCES` in `test_blastermaster.py`):

| File | Tolerance | Why |
|---|---|---|
| `vdw.vdw` | rtol = atol = 1e-4 | float32 rounding; values span 1e-13 to 1e24 |
| `trim.electrostatics.phi` | atol = 0.02 kT/e | iterative solver with platform math functions; values span ±1000 |
| `ligand.desolv.heavy`, `.hydrogen` | atol = 0.0011 | printed with 3 decimals: 1 in the last digit |

These cover the differences between the rebuilt and legacy programs, and between platforms
(which are much smaller). With the same compiler flags, different gcc versions produce identical
output, so on Linux the `rebuilt` controls are reproduced exactly.

For the legacy comparison, files that depend on sphgen's borderline decisions are skipped
(`SPHGEN_DERIVED_FILES`: its spheres, the spheres selected from them, and with thin spheres off
the electrostatics that use them). See
[03-cross-platform.md](03-cross-platform.md#differences-from-the-legacy-binaries) for the numbers.

## Updating the controls

When a change is *meant* to change blastermaster's output (a new reduce version, a bug fix in a
program, different defaults):

1. Make the change and run the tests: the failures show which files change, and by how much.
2. Look at the differences and make sure they are the intended ones.
3. Regenerate: `python tests/blastermaster/make_controls.py rebuilt` (runs both configs, ~20
   minutes). Use a Linux x86_64 build; the manylinux wheel and a local build with the same flags
   give identical results.
4. If `test_rebuilt_controls_match_legacy` now fails, decide whether the new differences from
   legacy are acceptable, and adjust `SPHGEN_DERIVED_FILES`/`TOLERANCES` only with a reason.
5. Commit the controls with a message that says why they changed.

Never regenerate controls just to make a failing test pass without understanding the difference.
The `legacy` controls should never be regenerated: they record what the original programs did.

## Other test fixtures

Only blastermaster is tested so far. dockopt, retrodock and the job schedulers have no tests
(`tests/` has only blastermaster's).
