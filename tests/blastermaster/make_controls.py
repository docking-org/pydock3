"""Generate control files that later blastermaster runs are compared against.

Usage: python make_controls.py <tier> [--out DIR]

`legacy` controls were generated with the original code and prebuilt binaries;
`rebuilt` controls with the binaries built from source by this package.
"""
import argparse
import shutil
import tempfile
from pathlib import Path

from harness import CONFIGS, CONTROLS_DIR, run_job


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tier", choices=["legacy", "rebuilt"])
    parser.add_argument("--out", type=Path, help="output dir (default: controls/<tier>)")
    args = parser.parse_args()
    out_dir = args.out or CONTROLS_DIR / args.tier

    for config_name, overrides in CONFIGS.items():
        config_out_dir = out_dir / config_name
        if config_out_dir.exists():
            shutil.rmtree(config_out_dir)
        config_out_dir.mkdir(parents=True)

        with tempfile.TemporaryDirectory() as tmp:
            produced = run_job(tmp, config_name)
            for name, path in sorted(produced.items()):
                shutil.copy(path, config_out_dir / name)
            shutil.copy(Path(tmp) / "job" / "blastermaster_config.yaml", config_out_dir)
        print(f"{config_name}: {len(produced)} files -> {config_out_dir}")


if __name__ == "__main__":
    main()
