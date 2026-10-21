"""Command-line manufactured solve and reproducible experiment suite."""

import argparse
import json
from pathlib import Path

from .experiments import run, write_results
from .fem import assemble, errors, solve, unit_square
from .problems import exact, gradient, source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    single = commands.add_parser("solve", help="solve the smooth manufactured problem")
    single.add_argument("--n", type=int, default=16)
    single.add_argument("--method", choices=["strong", "penalty", "nitsche"], default="nitsche")
    single.add_argument("--gamma", type=float, default=10)
    single.add_argument("--output", type=Path)
    study = commands.add_parser("study", help="reproduce CSV, JSON and figures")
    study.add_argument("--output", type=Path, default=Path("output/study"))
    study.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "study":
            write_results(args.output, run(), figures=not args.no_figures)
            print(f"Wrote measurements to {args.output.resolve()}")
        else:
            mesh = unit_square(args.n)
            solution = solve(assemble(mesh, source, exact), args.method, args.gamma)
            record = {
                "n": args.n,
                "method": args.method,
                "gamma": args.gamma,
                "gamma_star": solution.gamma_star,
                "backward_error": solution.backward_error,
                **errors(mesh, solution.values, exact, gradient),
            }
            text = json.dumps(record, indent=2, allow_nan=False) + "\n"
            if args.output:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(text, encoding="utf-8")
            print(text, end="")
    except (ValueError, ArithmeticError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
