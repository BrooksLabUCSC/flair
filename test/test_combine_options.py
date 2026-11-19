"""The input options of flair combine: --manifest names the isoform BEDs a line each,
--prefixes names the transcriptomes by prefix, --isoform_beds by path, and the two
latter take a list (cli_options)."""
import pytest
from flair import FlairInputDataError
from flair.flair_combine import parse_input


def write_beds(tmp_path, prefixes):
    for prefix in prefixes:
        bed = tmp_path / (prefix + '.isoforms.bed')
        bed.parent.mkdir(parents=True, exist_ok=True)
        bed.write_text('')


def test_prefixes_name_the_isoform_beds(tmp_path):
    write_beds(tmp_path, ('one', 'sub/two'))
    paths = parse_input(None, [f'{tmp_path}/one', f'{tmp_path}/sub/two'], None)
    assert paths == [f'{tmp_path}/one.isoforms.bed', f'{tmp_path}/sub/two.isoforms.bed']


def test_isoform_beds_are_paths():
    assert parse_input(None, None, ['a.bed,b.bed', 'c.bed']) == ['a.bed', 'b.bed', 'c.bed']


def test_manifest_names_one_bed_a_line(tmp_path):
    write_beds(tmp_path, ('one', 'two'))
    manifest = tmp_path / 'manifest.txt'
    manifest.write_text(f'# two samples\n{tmp_path}/one.isoforms.bed\n\n{tmp_path}/two.isoforms.bed\n')
    assert parse_input(str(manifest), None, None) == [f'{tmp_path}/one.isoforms.bed',
                                                      f'{tmp_path}/two.isoforms.bed']


def test_an_empty_manifest_is_an_error(tmp_path):
    manifest = tmp_path / 'manifest.txt'
    manifest.write_text('\n# nothing here\n')
    with pytest.raises(FlairInputDataError, match='no isoform BEDs'):
        parse_input(str(manifest), None, None)


def test_a_bed_named_twice_is_an_error(tmp_path):
    write_beds(tmp_path, ('one',))
    manifest = tmp_path / 'manifest.txt'
    manifest.write_text(f'{tmp_path}/one.isoforms.bed\n{tmp_path}/one.isoforms.bed\n')
    with pytest.raises(FlairInputDataError, match='named more than once'):
        parse_input(str(manifest), None, None)


def test_a_bed_the_manifest_names_must_exist(tmp_path):
    manifest = tmp_path / 'manifest.txt'
    manifest.write_text('nosuch.isoforms.bed\n')
    with pytest.raises(FlairInputDataError, match='isoform BED nosuch.isoforms.bed does not exist'):
        parse_input(str(manifest), None, None)
