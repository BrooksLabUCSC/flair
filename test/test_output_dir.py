"""The temporary directory a run makes for itself: kept only when it was asked for,
and removed whether the run finishes or fails."""
import os
import pytest
from flair import FlairInputDataError
from flair.io_utils import keep_or_remove, make_run_temp_dir, remove_on_failure


def test_a_run_directory_is_made_in_the_directory_given(tmp_path):
    run_dir = make_run_temp_dir(str(tmp_path / 'out'), str(tmp_path / 'tmp'))
    assert os.path.isdir(run_dir)
    assert run_dir.startswith(str(tmp_path / 'tmp'))
    assert 'out.' in os.path.basename(run_dir.rstrip('/'))


def test_a_run_directory_goes_to_the_system_one_when_none_is_given(tmp_path, monkeypatch):
    "tempfile.gettempdir() is what $TMPDIR sets, and what it caches for the process"
    import tempfile
    monkeypatch.setattr(tempfile, 'tempdir', str(tmp_path / 'systmp'))
    os.makedirs(tmp_path / 'systmp')
    run_dir = make_run_temp_dir(str(tmp_path / 'out'))
    assert os.path.isdir(run_dir) and run_dir.startswith(str(tmp_path / 'systmp'))


def test_a_finished_run_removes_its_directory(tmp_path):
    run_dir = make_run_temp_dir(str(tmp_path / 'out'), str(tmp_path / 'tmp'))
    with keep_or_remove(run_dir, False):
        pass
    assert not os.path.exists(run_dir)


def test_a_failed_run_removes_it_too(tmp_path):
    run_dir = make_run_temp_dir(str(tmp_path / 'out'), str(tmp_path / 'tmp'))
    with pytest.raises(FlairInputDataError):
        with keep_or_remove(run_dir, False):
            raise FlairInputDataError('nothing to test')
    assert not os.path.exists(run_dir)


def test_the_directory_asked_for_is_kept_and_reported(tmp_path, caplog):
    import logging
    run_dir = make_run_temp_dir(str(tmp_path / 'out'), str(tmp_path / 'tmp'))
    with caplog.at_level(logging.INFO):
        with keep_or_remove(run_dir, True):
            pass
    assert os.path.isdir(run_dir)
    assert 'intermediate files kept in' in caplog.text


def test_a_failed_run_leaves_no_directory_of_its_own(tmp_path):
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with pytest.raises(FlairInputDataError):
        with remove_on_failure(out_dir, True):
            raise FlairInputDataError('nothing to test')
    assert not os.path.exists(out_dir)


def test_a_failed_run_keeps_a_directory_it_was_given(tmp_path):
    "remove is False for a directory the run was given rather than made"
    out_dir = str(tmp_path / 'out')
    os.makedirs(out_dir)
    with pytest.raises(FlairInputDataError):
        with remove_on_failure(out_dir, False):
            raise FlairInputDataError('nothing to test')
    assert os.path.isdir(out_dir)
