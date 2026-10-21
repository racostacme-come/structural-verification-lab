"""Command-line analytical benchmark."""

import argparse
import json

from . import ContinuationError
from .validation import campaign


def main(argv=None):
    parser = argparse.ArgumentParser(description="Trace and validate a nonlinear two-bar truss")
    parser.add_argument(
        "--output", default="out", help="directory for CSV, JSON, and PNG outputs"
    )
    parser.add_argument(
        "--ratio", type=float, default=0.2, help="height/half-span in [0.01, 10]"
    )
    args = parser.parse_args(argv)
    try:
        summary = campaign(args.output, args.ratio)
    except (ValueError, ContinuationError, OSError, RuntimeError) as error:
        parser.exit(2, f"foldtrace: {error}\n")
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
