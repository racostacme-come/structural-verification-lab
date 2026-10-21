"""Command line entry point for a dimensionless column and sample study."""

import argparse
import json
from pathlib import Path

import numpy as np

from .model import solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("solve", help="unit-length, unit-EI column")
    run.add_argument("--elements", type=int, default=16)
    run.add_argument(
        "--support",
        choices=["pinned", "clamped", "cantilever", "clamped-pinned"],
        default="pinned",
    )
    run.add_argument(
        "--spring", type=float, default=0.0, help="equal end springs (pinned only)"
    )
    run.add_argument("--modes", type=int, default=1)
    run.add_argument("--output", type=Path)
    study = commands.add_parser("study", help="reproduce analytical validation and figures")
    study.add_argument("--output", type=Path, default=Path("output/study"))
    args = parser.parse_args()
    if args.command == "study":
        from .experiments import run_study

        run_study(args.output)
        return
    try:
        if args.elements < 1:
            raise ValueError("elements must be >= 1")
        result = solve(
            np.linspace(0, 1, args.elements + 1),
            support=args.support,
            springs=(args.spring, args.spring),
            modes=args.modes,
        )
    except (ValueError, ArithmeticError) as exc:
        parser.error(str(exc))
    payload = (
        json.dumps(
            {"loads": result.loads.tolist(), "residuals": result.residuals.tolist()}, indent=2
        )
        + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
