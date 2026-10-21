"""Command-line entry point."""

import argparse
import json
import sys

from . import __version__


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate and visualize adaptive 1D finite elements"
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--output", default="out", help="output directory (default: out)")
    args = parser.parse_args(argv)
    from .report import benchmark

    try:
        summary = benchmark(args.output)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"meshbeacon: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
