"""Compare blastermaster's output on the test fixtures with the control files."""
import pytest

from compare import Tolerance, compare_files
from harness import CONFIGS, CONTROLS_DIR, run_job, run_steps_in_isolation

TIER = "legacy"
TOLERANCES = {}  # file name -> Tolerance; exact by default
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


def control_file_names(control_dir):
    return {p.name for p in control_dir.iterdir()} - NOT_COMPARED


@pytest.mark.parametrize("config_name", CONFIGS)
def test_end_to_end(tmp_path, config_name):
    control_dir = CONTROLS_DIR / TIER / config_name
    produced = run_job(tmp_path, config_name)
    assert_matches_controls(produced, control_dir, control_file_names(control_dir))


@pytest.mark.parametrize("config_name", CONFIGS)
def test_steps_in_isolation(tmp_path, config_name):
    control_dir = CONTROLS_DIR / TIER / config_name
    produced = run_steps_in_isolation(tmp_path, config_name, control_dir)
    assert_matches_controls(produced, control_dir, control_file_names(control_dir) - {"INDOCK"})
