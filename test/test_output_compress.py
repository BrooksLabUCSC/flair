"""
The openers of flair's text outputs: the name a file is written under, and that a
file appears under that name only once it is complete.
"""
import pytest
from flair import FlairInputDataError
from flair.output_compress import (compressed_path, open_output, open_seekable_output,
                                   output_compression, set_output_compression)

OPENERS = (open_output, open_seekable_output)


@pytest.fixture
def compress(request):
    "a run writing its outputs compressed or plain, the setting restored after the test"
    was = output_compression()
    set_output_compression(request.param)
    yield request.param
    set_output_compression(was)


@pytest.mark.parametrize('compress', (False, True), indirect=True)
@pytest.mark.parametrize('opener', OPENERS)
def test_a_finished_output_is_installed_under_its_name(tmp_path, compress, opener):
    path = str(tmp_path / 'out.tsv')
    with opener(path) as fh:
        fh.write('a line\n')
    assert open(compressed_path(path), 'rb').read(1) == (b'\x1f' if compress else b'a')

@pytest.mark.parametrize('compress', (False, True), indirect=True)
@pytest.mark.parametrize('opener', OPENERS)
def test_a_failed_output_leaves_no_file_under_its_name(tmp_path, compress, opener):
    path = str(tmp_path / 'out.tsv')
    with pytest.raises(FlairInputDataError):
        with opener(path) as fh:
            fh.write('a line\n')
            raise FlairInputDataError('the run failed partway')
    assert list(tmp_path.iterdir()) == []
