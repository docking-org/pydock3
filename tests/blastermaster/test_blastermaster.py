"""Compare blastermaster's output on the test fixtures with the control files.

`rebuilt` controls: made with the programs built from source by this package (manylinux wheel).
`legacy` controls: made with the original programs (prebuilt with PGI compilers) and scripts.
"""
import pytest

from compare import Tolerance, compare_files
from harness import CONFIGS, CONTROLS_DIR, run_job, run_steps_in_isolation

# Floating-point differences between builds (file name -> tolerance; other files must match exactly)
TOLERANCES = {
    "vdw.vdw": Tolerance(rtol=1e-4, atol=1e-4),  # float32 rounding
    "trim.electrostatics.phi": Tolerance(atol=0.02),  # kT/e; values span about -1000 to 1000
    "ligand.desolv.heavy": Tolerance(atol=0.0011),  # 1 in the last printed digit
    "ligand.desolv.hydrogen": Tolerance(atol=0.0011),
}

# Files that depend on sphgen's borderline decisions, which differ between the legacy (PGI)
# build and ours for 0.5% of its spheres: its spheres and those selected from them, and the
# electrostatics computed with the low-dielectric spheres
SPHGEN_DERIVED_FILES = {
    "default": {"all_spheres.sph", "matching_spheres.sph", "matching_spheres.pdb"},
    "no_thin_spheres": {
        "all_spheres.sph", "matching_spheres.sph", "matching_spheres.pdb",
        "lowdielectric.sph", "lowdielectric.sph.pdb", "receptor.crg.lowdielectric.pdb",
        "qnifft.atm", "trim.electrostatics.phi",
    },
}

NOT_COMPARED = {"blastermaster_config.yaml"}


def assert_matches_controls(produced, control_dir, expected_file_names):
    missing = sorted(expected_file_names - produced.keys())
    unexpected = sorted(produced.keys() - expected_file_names)
    results = [
        compare_files(control_dir / name, produced[name], name, TOLERANCES.get(name, Tolerance()))
        for name in sorted(expected_file_names & produced.keys())
    ]
    problems = [str(r) for r in results if not r.ok]
    problems += [f"{name}: not produced" for name in missing]
    problems += [f"{name}: produced but not in controls" for name in unexpected]
    assert not problems, "\n".join(problems)


def control_files(control_dir):
    return {p.name: p for p in control_dir.iterdir() if p.name not in NOT_COMPARED}


@pytest.mark.parametrize("config_name", CONFIGS)
def test_end_to_end(tmp_path, config_name):
    control_dir = CONTROLS_DIR / "rebuilt" / config_name
    produced = run_job(tmp_path, config_name)
    assert_matches_controls(produced, control_dir, set(control_files(control_dir)))


@pytest.mark.parametrize("config_name", CONFIGS)
def test_steps_in_isolation(tmp_path, config_name):
    control_dir = CONTROLS_DIR / "rebuilt" / config_name
    produced = run_steps_in_isolation(tmp_path, config_name, control_dir)
    assert_matches_controls(produced, control_dir, set(control_files(control_dir)) - {"INDOCK"})


@pytest.mark.parametrize("config_name", CONFIGS)
def test_rebuilt_controls_match_legacy(config_name):
    rebuilt = control_files(CONTROLS_DIR / "rebuilt" / config_name)
    legacy_dir = CONTROLS_DIR / "legacy" / config_name
    comparable = set(control_files(legacy_dir)) - SPHGEN_DERIVED_FILES[config_name]
    assert_matches_controls(
        {name: path for name, path in rebuilt.items() if name in comparable}, legacy_dir, comparable
    )
