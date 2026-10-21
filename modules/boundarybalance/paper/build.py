"""Generate data-bound tables, compile LaTeX and enforce manuscript quality gates."""

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from pypdf import PdfReader


def sci(value):
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def numbers(root):
    data = json.loads((root / "examples/results/metrics.json").read_text())
    n32 = next(r for r in data["convergence"] if r["n"] == 32 and r["method"] == "nitsche")
    n8 = next(r for r in data["convergence"] if r["n"] == 8 and r["method"] == "nitsche")
    sweep = [r for r in data["sweep"] if r["method"] == "nitsche"]
    near = next(r for r in sweep if r["gamma"] == 1.57)
    high = next(r for r in sweep if r["gamma"] == 100000)
    commands = {
        "Critical": f"{data['gamma_star_n8']:.8f}",
        "OrderL": f"{data['orders']['nitsche']['l2']:.4f}",
        "OrderH": f"{data['orders']['nitsche']['h1']:.4f}",
        "PenaltyOrder": f"{data['orders']['penalty']['l2']:.4f}",
        "FineL": sci(n32["l2"]),
        "FineH": sci(n32["h1"]),
        "PatchN": sci(data["patch_max_error"]["nitsche"]),
        "PatchP": sci(data["patch_max_error"]["penalty"]),
        "Residual": sci(n8["backward_error"]),
        "CoarseL": sci(n8["l2"]),
        "NearL": sci(near["l2"]),
        "NearCondition": sci(near["condition_2"]),
        "HighCondition": sci(high["condition_2"]),
        "QuadLoad": sci(data["quadrature_check"]["load_max_difference"]),
        "QuadSolve": sci(data["quadrature_check"]["solution_max_difference"]),
        "QuadNorm": sci(data["quadrature_check"]["norm_max_difference"]),
    }
    rows = []
    for r in data["convergence"]:
        if r["method"] in ("nitsche", "penalty"):
            rows.append(
                f"{r['method'].capitalize()} & {r['n']} & ${sci(r['l2'])}$ & "
                f"${sci(r['h1'])}$ & ${sci(r['boundary_l2'])}$ " + r"\\"
            )
    commands["ConvergenceRows"] = "\n".join(rows)
    return (
        "\n".join(rf"\newcommand{{\{key}}}{{{value}}}" for key, value in commands.items())
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", help="path to tectonic or pdflatex")
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    root = folder.parent
    (folder / "numbers.tex").write_text(numbers(root), encoding="utf-8", newline="\n")
    engine = args.engine or shutil.which("tectonic") or shutil.which("pdflatex")
    if not engine:
        raise SystemExit("Install Tectonic or TeX Live; pass --engine if not on PATH")
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
        raise SystemExit("Supported engines: tectonic, pdflatex")
    for _ in range(passes):
        result = subprocess.run(
            command, cwd=folder, capture_output=True, text=True, errors="replace", timeout=180
        )
        if result.returncode:
            raise SystemExit(result.stdout + result.stderr)
    log = (build / "paper.log").read_text(errors="replace")
    issues = re.findall(
        r"Overfull \\[hv]box[^\n]*|[^\n]*(?:undefined references|"
        r"Citation .* undefined|Reference .* undefined)[^\n]*",
        log,
    )
    if issues:
        raise SystemExit("Layout/reference errors:\n" + "\n".join(issues))
    pdf = build / "paper.pdf"
    pages = PdfReader(pdf).pages
    if not 1 <= len(pages) <= 5:
        raise SystemExit(f"Page limit failed: {len(pages)}")
    if any(not (p.extract_text() or "").strip() for p in pages):
        raise SystemExit("Empty manuscript page")
    inputs = [
        folder / "paper.tex",
        folder / "numbers.tex",
        folder / "build.py",
        root / "pyproject.toml",
    ]
    inputs += sorted((root / "src").rglob("*.py"))
    inputs += sorted((root / "tests").glob("*.py"))
    inputs += sorted((root / "scripts").glob("*.py"))
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
