"""Mutation checks for published manuscript evidence and reproducibility."""
import json
import shutil

import pytest
from pypdf import PdfWriter
from manuscript import ROOT, FIXED, DATA, sources, verify, check_tables, build_report


@pytest.fixture
def manuscript_copy(tmp_path):
    for name in (*FIXED, *sources(ROOT), 'paper/paper.pdf', 'paper/build-report.json'):
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    return tmp_path


def test_published_manuscript_is_current():
    assert 1 <= verify() <= 5


@pytest.mark.parametrize('name', ['paper/paper.tex', 'paper/paper.pdf', 'suite.json'])
def test_modified_publication_is_rejected(manuscript_copy, name):
    with (manuscript_copy / name).open('ab') as stream:
        stream.write(b'\nchanged\n')
    with pytest.raises(ValueError, match='Stale manuscript evidence'):
        verify(manuscript_copy)


@pytest.mark.parametrize('change', ['remove_source', 'change_source', 'change_data', 'change_table'])
def test_stale_numerical_evidence_is_rejected(manuscript_copy, change):
    if change == 'remove_source':
        report_path = manuscript_copy / DATA[0]
        report = json.loads(report_path.read_text())
        report['source_sha256'].pop('native/assembly.cpp')
        report_path.write_text(json.dumps(report))
    elif change == 'change_source':
        with (manuscript_copy / 'native/assembly.cpp').open('ab') as stream:
            stream.write(b'\n// changed\n')
    elif change == 'change_data':
        report_path = manuscript_copy / DATA[0]
        report = json.loads(report_path.read_text())
        report['temporal_against_semidiscrete'][0]['force_error_N'] *= 2
        report_path.write_text(json.dumps(report))
    else:
        (manuscript_copy / 'paper/transfer-table.tex').write_text('stale table')
    with pytest.raises(ValueError, match='Stale'):
        check_tables(manuscript_copy)


@pytest.mark.parametrize('change', ['missing_hash', 'extra_hash', 'wrong_pages'])
def test_incomplete_manifest_is_rejected(manuscript_copy, change):
    path = manuscript_copy / 'paper/build-report.json'
    report = json.loads(path.read_text())
    if change == 'missing_hash':
        report['sha256'].pop('paper/paper.tex')
    elif change == 'extra_hash':
        report['sha256']['unused'] = '0' * 64
    else:
        report['pages'] = 0
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match='manifest|page count'):
        verify(manuscript_copy)


def test_six_page_manuscript_is_rejected(manuscript_copy):
    writer = PdfWriter()
    for _ in range(6):
        writer.add_blank_page(width=595, height=842)
    writer.write(manuscript_copy / 'paper/paper.pdf')
    with pytest.raises(ValueError, match='1-5 pages'):
        build_report(manuscript_copy)


def test_external_ci_pdf_and_manifest_are_bound(manuscript_copy, tmp_path):
    report = manuscript_copy / 'paper/build-report.json'
    pdf = manuscript_copy / 'paper/paper.pdf'
    assert verify(ROOT, report_path=report, pdf_path=pdf) == verify(ROOT)
    with pdf.open('ab') as stream:
        stream.write(b'\nchanged\n')
    with pytest.raises(ValueError, match='Stale manuscript evidence'):
        verify(ROOT, report_path=report, pdf_path=pdf)
