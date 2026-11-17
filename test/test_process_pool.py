"""
ProcessPool, the pool of worker processes the parallel subcommands run their tasks
in: what it does when a worker dies, when a task fails, and with the results.

Each task function is at module level so the forked workers can call it.
"""
import os
import signal
import pytest
from contextlib import contextmanager
from concurrent.futures.process import BrokenProcessPool
from flair import FlairInputDataError
from flair.process_pool import ProcessPool
from flair.pycbio.sys import fileOps

# longer than the tasks take, short enough to fail rather than stall a test run
HANG_TIMEOUT = 30


@contextmanager
def time_limit(seconds):
    "fail rather than wait longer than seconds, so a pool that hangs is a failed test"
    def on_alarm(signum, frame):
        raise AssertionError(f'the pool did not finish within {seconds} seconds')
    previous = signal.signal(signal.SIGALRM, on_alarm)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def _double(n):
    return 2 * n

def _kill_self_on(n):
    "SIGKILL in a worker, as the OOM killer sends it"
    if n == 2:
        os.kill(os.getpid(), signal.SIGKILL)
    return n

def _input_error_on(n):
    if n == 2:
        raise FlairInputDataError('a bad input')
    return n

def _record_run(packed):
    "write a file named for the task, so the tasks that ran can be counted"
    directory, n = packed
    if n == 0:
        raise FlairInputDataError('a bad input')
    open(os.path.join(directory, str(n)), 'w').close()
    return n


def test_imap_gives_the_results_in_the_order_of_the_tasks():
    with ProcessPool(3) as pool:
        assert list(pool.imap(_double, range(25))) == [2 * n for n in range(25)]

def test_run_runs_every_task(tmp_path):
    with ProcessPool(3) as pool:
        pool.run(_record_run, [(str(tmp_path), n) for n in range(1, 26)])
    assert sorted(int(p.name) for p in tmp_path.iterdir()) == list(range(1, 26))

def test_a_killed_worker_ends_the_run_rather_than_hanging():
    # multiprocessing.Pool waits forever for a worker that dies without a result
    with time_limit(HANG_TIMEOUT):
        with pytest.raises(BrokenProcessPool):
            with ProcessPool(3) as pool:
                pool.run(_kill_self_on, range(25))

def test_an_input_error_from_a_worker_keeps_no_worker_traceback():
    # the error reports without a stack, so the worker's traceback is not its cause
    with pytest.raises(FlairInputDataError) as info:
        with ProcessPool(3) as pool:
            pool.run(_input_error_on, range(25))
    assert info.value.__cause__ is None

def test_a_failed_task_leaves_most_of_the_tasks_unrun(tmp_path):
    # the first task fails, so only the ones already queued to the workers run
    with pytest.raises(FlairInputDataError):
        with ProcessPool(2) as pool:
            pool.run(_record_run, [(str(tmp_path), n) for n in range(500)])
    assert len(list(tmp_path.iterdir())) < 100

def test_a_file_opened_after_the_pool_is_made_is_not_held_by_a_worker(tmp_path):
    # a worker holding the pipe to bgzip open leaves it waiting for end of file
    path = str(tmp_path / 'out.tsv.gz')
    with time_limit(HANG_TIMEOUT):
        with ProcessPool(3) as pool:
            with fileOps.opengz(path, 'w', bgzip=True) as fh:
                fh.write('a line\n')
                assert list(pool.imap(_double, range(10))) == [2 * n for n in range(10)]
    assert os.path.getsize(path) > 0
