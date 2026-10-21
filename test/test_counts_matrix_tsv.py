"""
Tests for counts_matrix_tsv module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix_tsv import (CountsRow, ID_COLUMNS, CountsMatrixReader,
                                     read_counts_rows, read_sample_columns,
                                     read_isoform_ids, write_counts_matrix,
                                     describe_ids, MAX_REPORTED_IDS)

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
