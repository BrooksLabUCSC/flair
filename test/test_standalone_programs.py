"""The programs an install puts on PATH beside flair: what they report for a missing
input, and what they report about what they produced."""
import subprocess
import pytest

PROGRAMS = {
    'gtf_to_bed': (['gtf_to_bed', '/no/such.gtf', '{out}/o.bed'], 'annotation GTF file does not exist'),
    'bed_to_gtf': (['bed_to_gtf', '/no/such.bed'], 'isoform BED file does not exist'),
    'fasta_seq_lengths': (['fasta_seq_lengths', '/no/such.fa', '{out}/o.tsv'], 'FASTA file does not exist'),
    'mark_intron_retention': (['mark_intron_retention', '/no/such.bed', '{out}/o.bed', '{out}/o.txt'],
                              'isoform BED file does not exist'),
    'plot_isoform_usage': (['plot_isoform_usage', '/no/such.bed', '/no/such.tsv', 'KRAS'],
                           'isoform BED file does not exist'),
    'diff_iso_usage': (['diff_iso_usage', '/no/such.tsv', 'a', 'b', '{out}/o.tsv'],
                       'counts matrix file does not exist'),
    'diffsplice_fishers_exact': (['diffsplice_fishers_exact', '/no/such.tsv', 'a', 'b', '{out}/o.tsv'],
                                 'events quantification TSV file does not exist'),
    'spliceevents_outliers': (['spliceevents_outliers', '-i', '/no/such', '-o', '{out}/o'],
                              'splicing events counts TSV file does not exist'),
}


def run(cmd, tmp_path):
    return subprocess.run([arg.format(out=tmp_path) for arg in cmd],
                          capture_output=True, text=True)


@pytest.mark.parametrize('program', sorted(PROGRAMS))
def test_a_missing_input_is_named(program, tmp_path):
    "the message a subcommand's precheck gives, rather than the Python exception"
    cmd, expect = PROGRAMS[program]
    result = run(cmd, tmp_path)
    assert result.returncode == 1
    assert expect in result.stderr
    assert len(result.stderr.strip().split('\n')) == 2, result.stderr


def test_a_gtf_with_no_transcripts_is_a_warning(tmp_path):
    "a BED given where a GTF was meant wrote an empty BED and reported nothing"
    gtf = tmp_path / 'not.gtf'
    gtf.write_text('chr1\t100\t200\tiso1\t0\t+\n')
    result = run(['gtf_to_bed', str(gtf), '{out}/o.bed'], tmp_path)
    assert result.returncode == 0
    assert 'no transcripts' in result.stderr


def test_a_fasta_with_no_sequences_is_a_warning(tmp_path):
    fasta = tmp_path / 'not.fa'
    fasta.write_text('nothing here\n')
    result = run(['fasta_seq_lengths', str(fasta), '{out}/o.tsv'], tmp_path)
    assert result.returncode == 0
    assert 'no sequences' in result.stderr


def test_what_was_written_is_reported(tmp_path):
    gtf = tmp_path / 'one.gtf'
    gtf.write_text('chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    result = run(['gtf_to_bed', str(gtf), '{out}/o.bed'], tmp_path)
    assert result.returncode == 0
    assert '1 transcript' in result.stderr
