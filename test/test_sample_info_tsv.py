"""
Tests for sample_info_tsv module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix_tsv import CountsRow, write_counts_matrix
from flair.sample_info_tsv import (SampleInfo, sample_info_path, write_sample_info,
                                   read_sample_info)

SAMPLE_INFOS = (SampleInfo('s1', 'ctl', 'b1'),
                SampleInfo('s2', 'ctl', 'b1'),
                SampleInfo('s3', 'test', 'b2'))

def write_test_matrix(path, sample_columns):
    write_counts_matrix(path, sample_columns,
                        [CountsRow('gene1', 'iso1', [1] * len(sample_columns))])

def test_sample_info_path_counts_tsv():
    assert sample_info_path('out/flair.quantify.counts.tsv') == 'out/flair.quantify.sample_info.tsv'

def test_sample_info_path_other_names():
    assert sample_info_path('matrix.tsv') == 'matrix.sample_info.tsv'
    assert sample_info_path('matrix') == 'matrix.sample_info.tsv'

def test_read_sample_info_round_trip(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, [si.sample_id for si in SAMPLE_INFOS])
    write_sample_info(sample_info_path(counts), SAMPLE_INFOS)
    assert read_sample_info(counts) == list(SAMPLE_INFOS)

def test_read_sample_info_allows_underscores(tmp_path):
    "the fields are columns of their own, so they may contain the joining character"
    infos = [SampleInfo('s_1', 'ctl_a', 'b_1'), SampleInfo('s_2', 'test_a', 'b_2')]
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, [si.sample_id for si in infos])
    write_sample_info(sample_info_path(counts), infos)
    assert read_sample_info(counts) == infos

def test_read_sample_info_falls_back_to_column_names(tmp_path):
    "a counts matrix written before the sample info file existed"
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, ['s1_ctl_b1', 's2_test_b1'])
    assert read_sample_info(counts) == [SampleInfo('s1_ctl_b1', 'ctl', 'b1'),
                                        SampleInfo('s2_test_b1', 'test', 'b1')]

def test_read_sample_info_must_match_counts_columns(tmp_path):
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, ['s1', 's2'])
    write_sample_info(sample_info_path(counts), SAMPLE_INFOS)
    with pytest.raises(FlairInputDataError, match="does not describe the columns"):
        read_sample_info(counts)

def test_sample_info_wrong_columns(tmp_path):
    sample_info = str(tmp_path / 'x.sample_info.tsv')
    with open(sample_info, 'w') as fh:
        fh.write('sample\tcondition\tbatch\ns1\tctl\tb1\n')
    counts = str(tmp_path / 'x.counts.tsv')
    write_test_matrix(counts, ['s1'])
    with pytest.raises(FlairInputDataError, match="expected columns"):
        read_sample_info(counts)
