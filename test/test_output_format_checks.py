"""The checks that keep what a run writes readable: the sample names that go into
the files, the ids of a counts matrix, and a read map of a FLAIR before 3.1."""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix_tsv import CountsRow, check_isoform_ids_unique, check_sample_columns
from flair.read_map_tsv import ReadMapReader, ReadMapWriter
from flair.read_map_index import ReadMapIndex
from flair.run_checks import check_sample_name


@pytest.mark.parametrize('name, what', [('a,b', 'comma'), ('a\tb', 'tab'), ('a/b', 'slash')])
def test_a_sample_name_that_the_files_would_not_survive(name, what):
    with pytest.raises(FlairInputDataError, match=what):
        check_sample_name(name)


def test_a_sample_name_is_otherwise_the_user_s_own():
    for name in ('sample 1', 'A1', 'tumor.rep-2', "o'brien"):
        check_sample_name(name)


def test_the_manifest_names_the_file_the_name_came_from():
    with pytest.raises(FlairInputDataError, match='manifest.tsv: sample'):
        check_sample_name('a,b', where='manifest.tsv: sample')


def test_a_repeated_isoform_id_is_reported():
    rows = [CountsRow('G1', 'I1', [1]), CountsRow('G1', 'I2', [1]), CountsRow('G1', 'I1', [2])]
    with pytest.raises(FlairInputDataError, match='1 repeated isoform_id: I1'):
        check_isoform_ids_unique('counts.tsv', rows)


def test_ids_that_appear_once_are_taken():
    check_isoform_ids_unique('counts.tsv', [CountsRow('G1', 'I1', [1]), CountsRow('G1', 'I2', [1])])


def test_two_columns_of_one_sample_are_reported():
    with pytest.raises(FlairInputDataError, match='sample A1 names more than one column'):
        check_sample_columns('counts.tsv', ['A1', 'B1', 'A1'])


def test_a_read_map_round_trips_under_its_column_name(tmp_path):
    path = str(tmp_path / 'x.read.map.tsv')
    with ReadMapWriter(path) as writer:
        writer.writeReads('FLT00000001|H1', ['r1', 'r2'])
    assert open(path).read().split('\n')[0] == 'isoform_id\treads'
    assert [(row.isoform_id, row.reads) for row in ReadMapReader(path)] == [('FLT00000001|H1', ['r1', 'r2'])]


def _old_read_map(tmp_path):
    path = str(tmp_path / 'old.read.map.tsv')
    with open(path, 'w') as fh:
        fh.write('name\treads\nFLT00000001\tr1,r2\n')
    return path


def test_the_read_map_of_an_earlier_flair_is_reported(tmp_path):
    with pytest.raises(FlairInputDataError, match='a FLAIR before 3.1'):
        ReadMapReader(_old_read_map(tmp_path))


def test_the_index_reports_it_too(tmp_path):
    "indexed by offset, so its first line would otherwise be an isoform called name"
    with pytest.raises(FlairInputDataError, match='a FLAIR before 3.1'):
        ReadMapIndex(_old_read_map(tmp_path), str(tmp_path))
