"""The output directory of flair diffexp and flair diffsplice: made by the run, kept
only when the run wrote something, and never shared with another run's files."""
import os
import pytest
from flair import FlairInputDataError
from flair.io_utils import remove_on_failure
from flair.flair_diffExp import make_output_dir as diffexp_output_dir
from flair.flair_diffSplice import make_output_dir as diffsplice_output_dir

MAKE_OUTPUT_DIRS = (diffexp_output_dir, diffsplice_output_dir)


@pytest.mark.parametrize('make_output_dir', MAKE_OUTPUT_DIRS)
def test_a_new_directory_is_made_with_its_workdir(tmp_path, make_output_dir):
    out_dir = str(tmp_path / 'out')
    workdir = os.path.join(out_dir, 'workdir')
    assert make_output_dir(out_dir, workdir, False)
    assert os.path.isdir(workdir)


@pytest.mark.parametrize('make_output_dir', MAKE_OUTPUT_DIRS)
def test_a_directory_that_is_there_is_an_error(tmp_path, make_output_dir):
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with pytest.raises(FlairInputDataError, match='already exists'):
        make_output_dir(out_dir, os.path.join(out_dir, 'workdir'), False)


@pytest.mark.parametrize('make_output_dir', MAKE_OUTPUT_DIRS)
def test_overwrite_takes_a_directory_that_is_there(tmp_path, make_output_dir):
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    assert not make_output_dir(out_dir, os.path.join(out_dir, 'workdir'), True)
    assert os.path.isdir(os.path.join(out_dir, 'workdir'))


def test_a_failed_run_leaves_no_directory_of_its_own(tmp_path):
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with pytest.raises(FlairInputDataError):
        with remove_on_failure(out_dir, True):
            raise FlairInputDataError('nothing to test')
    assert not os.path.exists(out_dir)


def test_a_failed_run_keeps_a_directory_it_was_given(tmp_path):
    "with --overwrite_output the directory is the user's, holding another run's files"
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with pytest.raises(FlairInputDataError):
        with remove_on_failure(out_dir, False):
            raise FlairInputDataError('nothing to test')
    assert os.path.isdir(out_dir)


def test_a_run_that_finishes_keeps_its_directory(tmp_path):
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with remove_on_failure(out_dir, True):
        pass
    assert os.path.isdir(out_dir)
