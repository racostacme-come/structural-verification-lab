"""Generate measured tables, compile LaTeX, and reject oversized or broken manuscripts."""

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from pypdf import PdfReader


def scientific(value):
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", help="Tectonic or pdfLaTeX executable")
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    root = folder.parent
    data = json.loads((root / "examples/results/metrics.json").read_text())
    commands = []
    for key, macro in [("0.0", "PinnedOrder"), ("10.0", "SpringOrder")]:
        commands.append(rf"\newcommand{{\{macro}}}{{{data['final_orders'][key]:.4f}}}")
    spring = next(r for r in data["restraint"] if r["kappa"] == 10)
    commands += [
        rf"\newcommand{{\SpringLoad}}{{{spring['load']:.8f}}}",
        rf"\newcommand{{\SpringGain}}{{{spring['exact'] / data['boundary'][0]['exact']:.5f}}}",
    ]
    coarse = next(r for r in data["convergence"] if r["kappa"] == 10 and r["elements"] == 4)
    commands += [
        rf"\newcommand{{\CoarsePercent}}{{{100 * coarse['relative_load_error']:.4f}}}",
        rf"\newcommand{{\CoarseResidual}}{{{scientific(coarse['backward_error'])}}}",
    ]
    rows = []
    for r in data["boundary"]:
        rows.append(
            f"{r['support']} & {r['exact']:.8f} & {r['load']:.8f} & "
            f"${scientific(r['relative_error'])}$ " + r"\\"
        )
    commands.append(r"\newcommand{\BoundaryRows}{" + "\n".join(rows) + "}")
    rows = []
    for r in data["convergence"]:
        rows.append(
            f"{r['kappa']:g} & {r['elements']} & "
            f"${scientific(r['relative_load_error'])}$ & "
            f"${scientific(r['field_l2_error'])}$ & "
            f"${scientific(r['backward_error'])}$ " + r"\\"
        )
    commands.append(r"\newcommand{\ConvergenceRows}{" + "\n".join(rows) + "}")
    (folder / "numbers.tex").write_text("\n".join(commands) + "\n", newline="\n")
    engine = args.engine or shutil.which("tectonic") or shutil.which("pdflatex")
    if not engine:
        raise SystemExit("Install Tectonic or TeX Live, then pass --engine.")
    if Path(engine).exists():
        engine = str(Path(engine).resolve())
    name = Path(engine).stem.lower()
    build = folder / "build"
    build.mkdir(exist_ok=True)
    if name == "tectonic":
        command = [
            engine,
            "--keep-logs",
            "--keep-intermediates",
            "--outdir",
            str(build),
            "paper.tex",
        ]
        passes = 1
    elif name == "pdflatex":
        command = [
            engine,
            "-no-shell-escape",
            "-halt-on-error",
            "-interaction=nonstopmode",
            "-output-directory",
            str(build),
            "paper.tex",
        ]
        passes = 2
    else:
        raise SystemExit("Only tectonic and pdflatex are supported.")
    for _ in range(passes):
        result = subprocess.run(
            command, cwd=folder, capture_output=True, text=True, errors="replace", timeout=180
        )
        if result.returncode:
            raise SystemExit(result.stdout + result.stderr)
    log = (build / "paper.log").read_text(errors="replace")
    failures = re.findall(
        r"Overfull \\[hv]box[^\n]*|[^\n]*(?:undefined references|"
        r"Citation .* undefined|Reference .* undefined)[^\n]*",
        log,
    )
    if failures:
        raise SystemExit("Layout/reference errors:\n" + "\n".join(failures))
    pdf = build / "paper.pdf"
    pages = PdfReader(pdf).pages
    if not 1 <= len(pages) <= 5:
        raise SystemExit(f"Page limit failed: {len(pages)} pages")
    if any(not (page.extract_text() or "").strip() for page in pages):
        raise SystemExit("Empty manuscript page")
    inputs = [
        folder / "paper.tex",
        folder / "numbers.tex",
        folder / "build.py",
        root / "pyproject.toml",
        root / "CMakeLists.txt",
    ]
    inputs += sorted((root / "src").rglob("*.py"))
    inputs += sorted((root / "cpp").glob("*"))
    inputs += sorted((root / "examples/results").glob("*"))
    report = {
        "pages": len(pages),
        "page_limit": 5,
        "engine": name,
        "overfull_boxes": 0,
        "undefined_references": 0,
        "inputs_sha256": {
            p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs
        },
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
    }
    shutil.copyfile(pdf, folder / "paper.pdf")
    (folder / "validation.json").write_text(json.dumps(report, indent=2) + "\n", newline="\n")
    print(f"Compiled {len(pages)} pages; no overfull boxes or undefined references.")


if __name__ == "__main__":
    main()
