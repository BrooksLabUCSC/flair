"""
ReadMapIndex: an isoform's reads from a read map, plain, bgzip or gzip, without
reading the rest
"""
import gzip
import subprocess
import pytest
from flair import FlairInputDataError
from flair.read_map_index import ReadMapIndex, is_bgzip

ROWS = 'name\treads\niso1\tr1,r2\niso2\tr3\niso3\tr4,r5,r6\n'


def _write(path, text=ROWS):
    with open(path, 'w') as fh:
        fh.write(text)
    return str(path)


def _assigned(path, tmp_path, isoforms):
    with ReadMapIndex(path, str(tmp_path)).open_reader() as reader:
        return reader.assigned(isoforms)


@pytest.mark.parametrize('compress', ['plain', 'bgzip', 'gzip'])
def test_only_the_isoforms_asked_for_are_read(tmp_path, compress):
    path = _write(tmp_path / 'x.read.map.tsv')
    if compress == 'bgzip':
        subprocess.run(['bgzip', path], check=True)
        path += '.gz'
        assert is_bgzip(path)
    elif compress == 'gzip':
        with open(path, 'rb') as src, gzip.open(path + '.gz', 'wb') as dst:
            dst.write(src.read())
        path += '.gz'
        assert not is_bgzip(path)
    assert _assigned(path, tmp_path, ['iso3', 'iso1', 'missing']) == {'r4': 'iso3', 'r5': 'iso3', 'r6': 'iso3',
                                                                      'r1': 'iso1', 'r2': 'iso1'}


def test_a_map_without_a_header_is_read_too(tmp_path):
    path = _write(tmp_path / 'x.read.map.tsv', ROWS.split('\n', 1)[1])
    assert _assigned(path, tmp_path, ['iso2']) == {'r3': 'iso2'}


def test_an_isoform_of_two_rows_is_an_error(tmp_path):
    path = _write(tmp_path / 'x.read.map.tsv', 'iso1\tr1\niso1\tr2\n')
    with pytest.raises(FlairInputDataError, match='two rows'):
        ReadMapIndex(path, str(tmp_path))
