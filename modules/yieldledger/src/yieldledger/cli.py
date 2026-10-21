"""Command-line entry point for the reproducible validation suite."""

import argparse
import json

from .experiments import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results", help="directory for generated artifacts")
    args = parser.parse_args()
    print(json.dumps(run(args.output), indent=2))


if __name__ == "__main__":
    main()
