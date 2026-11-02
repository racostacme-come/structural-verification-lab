"""Bind the thermoelastic manuscript to measured JSON and versioned source.

Hashes establish identity, not numerical or physical validity. Reproduction
and independent numerical gates remain in tests/test_* and the CI suite job.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ('results/coupled-refinement.json', 'results/coupled-transfer.json')
TABLES = ('paper/refinement-table.tex', 'paper/transfer-table.tex',
          'paper/contact-table.tex', 'paper/measurements.tex')
FIXED = ('paper/paper.tex', 'manuscript.py', 'suite.py', 'suite.json', *DATA, *TABLES)


def sources(root):
    return sorted(p.relative_to(root).as_posix()
                  for folder, pattern in (('thermomechanics', '*.py'), ('native', '*'))
                  for p in (root / folder).glob(pattern) if p.is_file())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_evidence(root=ROOT):
    root = Path(root)
    expected = {name: hashlib.sha256((root / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                for name in sources(root)}
    if not expected:
        raise ValueError('Missing numerical sources')
    data = [json.loads((root / name).read_text(encoding='utf-8')) for name in DATA]
    for name, report in zip(DATA, data, strict=True):
        if report.get('source_sha256') != expected:
            raise ValueError(f'Stale or incomplete numerical source manifest: {name}')
    return data


def table(columns, header, rows):
    return '\n'.join([r'\begin{center}\small', r'\begin{tabular}{' + columns + r'}\toprule',
                      header + r'\\\midrule', *[r + r'\\' for r in rows],
                      r'\bottomrule\end{tabular}', r'\end{center}', ''])


def render_tables(refinement, coupled):
    rows = []
    for label, key, count in (
        ('Time', 'temporal_against_semidiscrete', 'steps'),
        ('Space', 'spatial_exact_semidiscrete_against_continuum', 'cells'),
    ):
        for row in refinement[key]:
            order = row.get('orders', {}).get('temperature_error_K')
            rows.append(f"{label} & {row[count]} & {row['temperature_error_K']:.3e} & "
                        f"{row['force_error_N']:.3e} & {row['nodal_error_m']:.3e} & "
                        + ('--' if order is None else f'{order:.4f}'))
    refinement_table = table('lrrrrr',
        r'Study & Steps/cells & $e_T$ (K) & $e_N$ (N) & $e_u$ (m) & $p_T$', rows)
    transfer_table = table('rrrr',
        r'Elements & Integral error & Overlap force error & Naive force error',
        [f"{r['elements']} & {r['scaled_integral_error']:.3e} & "
         f"{r['overlap_relative_force_error']:.3e} & {r['naive_relative_force_error']:.3e}"
         for r in coupled['transfer_only']])
    contact = refinement['contact_property_sweep']
    contact_rows = []
    for resistance in sorted({r['resistance_m2_K_per_W'] for r in contact}):
        group = [r for r in contact if r['resistance_m2_K_per_W'] == resistance]
        contact_rows.append(f'{resistance:.1e} & ' + ' & '.join(
            f'{max(r[key] for r in group):.3e}' for key in
            ('relative_rate_error', 'relative_force_error', 'nodal_error_m')))
    contact_table = table('rrrr',
        r'$R$ (m$^2$ K/W) & Heat-rate rel. error & Force rel. error & $e_u$ (m)', contact_rows)
    stats = {
        'MaxIntegralError': max(r['scaled_integral_error'] for r in coupled['transfer_only']),
        'MaxEnergyError': coupled['maximum_scaled_errors']['scaled_energy_residual'],
        'MaxSourceForceError': coupled['maximum_scaled_errors']['scaled_force_error'],
        'MaxContactDisplacement': max(r['nodal_error_m'] for r in contact),
        'MaxContactCondition': max(r['heat_condition_number'] for r in contact),
    }
    macros = '\n'.join('\\newcommand{\\' + key + '}{' + f'{value:.4e}' + '}'
                       for key, value in stats.items()) + '\n'
    return dict(zip(TABLES, (refinement_table, transfer_table, contact_table, macros), strict=True))


def check_tables(root=ROOT, write=False):
    root = Path(root)
    for name, expected in render_tables(*read_evidence(root)).items():
        path = root / name
        if write:
            path.write_text(expected, encoding='utf-8', newline='\n')
        elif not path.is_file() or path.read_text(encoding='utf-8') != expected:
            raise ValueError(f'Stale manuscript table: {name}; run python manuscript.py --write-tables')


def build_report(root=ROOT):
    from pypdf import PdfReader
    root = Path(root)
    check_tables(root)
    pages = len(PdfReader(root / 'paper/paper.pdf').pages)
    if not 1 <= pages <= 5:
        raise ValueError(f'Manuscript must have 1-5 pages, got {pages}')
    names = (*FIXED, *sources(root), 'paper/paper.pdf')
    return {'pages': pages, 'sha256': {name: digest(root / name) for name in names}}


def verify(root=ROOT, report_path=None, pdf_path=None):
    from pypdf import PdfReader
    root = Path(root)
    check_tables(root)
    report_path = Path(report_path) if report_path else root / 'paper/build-report.json'
    pdf_path = Path(pdf_path) if pdf_path else root / 'paper/paper.pdf'
    report = json.loads(report_path.read_text(encoding='utf-8'))
    names = {*FIXED, *sources(root), 'paper/paper.pdf'}
    if set(report.get('sha256', {})) != names:
        raise ValueError('Manuscript manifest must include exactly every required input and PDF')
    for name in names:
        if digest(pdf_path if name == 'paper/paper.pdf' else root / name) != report['sha256'][name]:
            raise ValueError(f'Stale manuscript evidence: {name}')
    pages = len(PdfReader(pdf_path).pages)
    if not 1 <= pages <= 5 or report.get('pages') != pages:
        raise ValueError(f'Manuscript page count mismatch or outside 1-5: {pages}')
    return pages


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-tables', action='store_true')
    parser.add_argument('--report', type=Path)
    parser.add_argument('--pdf', type=Path)
    args = parser.parse_args()
    if args.write_tables:
        check_tables(write=True)
        print('Generated four manuscript evidence files from source-checked JSON.')
    else:
        print(f'All manuscript hashes and {verify(report_path=args.report, pdf_path=args.pdf)} pages verified.')
