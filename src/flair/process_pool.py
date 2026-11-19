"""Running tasks in a pool of worker processes.

The workers are forked, so each inherits what the parent built, such as the
annotation of flair transcriptome or the lookups of flair quantify, rather than
being sent it.

A ProcessPoolExecutor is used rather than multiprocessing.Pool, which waits
forever for a worker that dies without returning a result, as a worker killed by
the OOM killer does; the executor raises BrokenProcessPool instead.

Tasks are submitted a few per worker ahead of the workers rather than all at
once, since the executor has no way to stop a task it has already queued to
them.  After a task fails, the tasks running and the one or two rounds the
executor has queued still run, rather than all of the tasks left.  Leaving the
pool waits for those, rather than killing the workers and leaving the programs
they are running, such as minimap2, running after the run has gone.
"""
import logging
import multiprocessing as mp
from collections import deque
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from flair.pycbio import NoStackError

# tasks queued per worker, enough that a worker never waits for the main process
QUEUED_PER_WORKER = 4


def _fork_workers():
    "the task the workers are forked with, which does nothing else"


def fork_context():
    """The multiprocessing context of the pools and of the state shared with them.
    Fork, so a worker inherits what the parent built: pickling an annotation into
    every task took longer than the tasks.  Python 3.14 makes forkserver the default
    on Linux, so this is specified rather than assumed."""
    return mp.get_context('fork')


def _task_result(future):
    """The result of one task.  An error from a worker comes back with the worker's
    traceback attached as its cause, which for an error that reports without a stack,
    an input error, is the detail the report leaves out, so it goes to the debug log
    and the error is raised with the cause dropped."""
    try:
        return future.result()
    except Exception as exc:
        if not isinstance(exc, NoStackError):
            raise
        logging.debug('error in a pool worker', exc_info=True)
        raise exc from None


def _take_finished(pending):
    "wait for tasks to finish, then take their results and drop them from pending"
    done, _ = wait(pending, return_when=FIRST_COMPLETED)
    for future in done:
        pending.remove(future)
        _task_result(future)


class ProcessPool:
    """A pool of worker processes, threads of them, running the tasks passed to run()
    or imap().

    Used in a with statement, which on the way out drops the tasks not yet queued to a
    worker and waits for the rest.
    """
    def __init__(self, threads, *, initializer=None, initargs=()):
        self._queued = QUEUED_PER_WORKER * threads
        self._executor = ProcessPoolExecutor(threads, mp_context=fork_context(),
                                             initializer=initializer, initargs=initargs)
        # the workers are forked here, with a task that does nothing, rather than at
        # the first task: a worker that inherits a file the caller opens after the
        # pool is made keeps it open, and a pipe to a compressing program then never
        # reaches end of file
        self._executor.submit(_fork_workers).result()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._executor.shutdown(cancel_futures=True)
        return False

    def imap(self, func, tasks):
        "the result of func for each task, in the order of the tasks"
        pending = deque()
        for task in tasks:
            pending.append(self._executor.submit(func, task))
            if len(pending) > self._queued:
                yield _task_result(pending.popleft())
        while pending:
            yield _task_result(pending.popleft())

    def run(self, func, tasks):
        """Run func on each task, discarding the results.  A result is taken as its
        task finishes rather than in the order of the tasks, so a long task does not
        block the queue while the other workers idle."""
        pending = []
        for task in tasks:
            pending.append(self._executor.submit(func, task))
            if len(pending) > self._queued:
                _take_finished(pending)
        while pending:
            _take_finished(pending)
