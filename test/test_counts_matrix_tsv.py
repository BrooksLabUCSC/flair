"""
Tests for counts_matrix_tsv module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix_tsv import (CountsRow, ID_COLUMNS, CountsMatrixReader,
                                     read_counts_rows, read_sample_columns,
                                     read_isoform_ids, write_counts_matrix,
                                     describe_ids, MAX_REPORTED_IDS, read_grouped_counts,
                                     QUANTIFY_COUNTS, SPLICEEVENTS_COUNTS, VARIANTQUANT_COUNTS)
from flair import splicing_events_tsv, var_counts_tsv

ROWS = (CountsRow('gene1', 'iso1', [1, 2]), CountsRow('gene1', 'iso2', [3, 4]))

def write_old_format(path):
    with open(path, 'w') as fh:
        fh.write('ids\ts1\ts2\n')
        fh.write('iso1_gene1\t1\t2\n')

def test_round_trip(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_counts_matrix(counts, ['s1', 's2'], ROWS)
    assert read_sample_columns(counts) == ['s1', 's2']
    assert read_counts_rows(counts) == list(ROWS)
    assert read_isoform_ids(counts) == ['iso1', 'iso2']

def test_round_trip_compressed(tmp_path):
    "the suffix chooses the compression, nothing else changes"
    counts = str(tmp_path / 'x.counts.tsv.gz')
    write_counts_matrix(counts, ['s1', 's2'], ROWS)
    assert open(counts, 'rb').read(2) == b'\x1f\x8b'
    assert read_counts_rows(counts) == list(ROWS)

def test_counts_are_ints(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_counts_matrix(counts, ['s1', 's2'], ROWS)
    assert read_counts_rows(counts)[0].counts == [1, 2]

def test_non_numeric_count_is_an_error(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    with open(counts, 'w') as fh:
        fh.write('\t'.join(list(ID_COLUMNS) + ['s1']) + '\n')
        fh.write('gene1\tiso1\tnope\n')
    with pytest.raises(Exception, match="Error reading TSV row"):
        read_counts_rows(counts)

def test_old_format_refused(tmp_path):
    "a matrix with one isoform_gene column cannot be split reliably, so it is refused"
    counts = str(tmp_path / 'x.counts.tsv')
    write_old_format(counts)
    with pytest.raises(FlairInputDataError, match="must start with the columns"):
        CountsMatrixReader(counts)

def test_wrong_count_of_counts(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    with pytest.raises(FlairInputDataError, match="has 1 counts"):
        write_counts_matrix(counts, ['s1', 's2'], [CountsRow('gene1', 'iso1', [1])])

def test_failed_write_leaves_no_file(tmp_path):
    "the writer is atomic, so a crash part way through installs nothing"
    counts = str(tmp_path / 'x.counts.tsv')
    with pytest.raises(FlairInputDataError):
        write_counts_matrix(counts, ['s1'], [CountsRow('gene1', 'iso1', [1, 2])])
    assert not (tmp_path / 'x.counts.tsv').exists()

def test_describe_ids_elides():
    ids = ['id%d' % num for num in range(MAX_REPORTED_IDS + 1)]
    assert describe_ids(ids).endswith(', ...')
    assert describe_ids(ids[:2]) == 'id0, id1'


# read_grouped_counts: a quantify matrix, or the counts of spliceevents or variantquant
def _write_spliceevents_counts(path, rows):
    "rows (event, alternative, counts...) of a spliceevents counts table of samples s1, s2"
    with open(path, 'w') as fh:
        fh.write('\t'.join(list(splicing_events_tsv.EVENT_COLUMNS) + ['s1', 's2']) + '\n')
        for event, alternative, *counts in rows:
            fh.write('\t'.join([event, alternative, 'es', 'G1'] + [''] * 8 + [str(c) for c in counts]) + '\n')


def test_a_quantify_matrix_keeps_its_ids(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_counts_matrix(counts, ['s1', 's2'], ROWS)
    grouped = read_grouped_counts(counts)
    assert (grouped.format, grouped.sample_columns, grouped.rows, grouped.unreported) == \
        (QUANTIFY_COUNTS, ['s1', 's2'], list(ROWS), set())


def test_spliceevents_rows_are_grouped_by_event_and_other_is_unreported(tmp_path):
    counts = str(tmp_path / 'x.spliceevents.counts.tsv')
    event = 'es|chr1:101-200(+)|G1'
    _write_spliceevents_counts(counts, [(event, 'inc', 5, 1), (event, 'exc', 2, 6), (event, 'other', 1, 0)])
    grouped = read_grouped_counts(counts)
    assert grouped.format is SPLICEEVENTS_COUNTS
    assert read_sample_columns(counts) == grouped.sample_columns == ['s1', 's2']
    assert grouped.rows == [CountsRow(event, event + '|inc', [5, 1]), CountsRow(event, event + '|exc', [2, 6]),
                            CountsRow(event, event + '|other', [1, 0])]
    assert grouped.unreported == {event + '|other'}


def test_variantquant_rows_are_grouped_by_var_name(tmp_path):
    counts = str(tmp_path / 'x.var.counts.tsv')
    with var_counts_tsv.var_counts_writer(counts, ['s1', 's2']) as writer:
        writer.writeAllele('chr1:5', 'A', 'G', ('.', '100', '2'), 'G1', [(7, 3), (4, 0)])
    grouped = read_grouped_counts(counts)
    assert grouped.format is VARIANTQUANT_COUNTS
    assert grouped.rows == [CountsRow('chr1:5:A>G|100|2', 'chr1:5:A>G|100|2|yes', [3, 0]),
                            CountsRow('chr1:5:A>G|100|2', 'chr1:5:A>G|100|2|no', [7, 4])]


def test_a_table_of_psi_values_is_refused(tmp_path):
    psi = str(tmp_path / 'x.spliceevents.PSIjunc.tsv')
    _write_spliceevents_counts(psi, [('es|chr1:101-200(+)|G1', 'inc', 0.5, '')])
    with pytest.raises(FlairInputDataError, match='must hold read counts'):
        read_grouped_counts(psi)


def test_an_unknown_table_is_refused(tmp_path):
    other = str(tmp_path / 'x.tsv')
    with open(other, 'w') as fh:
        fh.write('name\tkind\ts1\nA\tB\t1\n')
    with pytest.raises(FlairInputDataError, match='not a table of counts'):
        read_grouped_counts(other)
