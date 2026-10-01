"""Share a strict thread budget between parallel tasks and the programs they run.

Tasks such as transcriptome partitions or quantify genes run in a
multiprocessing pool with one worker per thread, and each running task counts
as one thread.  A program a task runs, such as minimap2, gets one thread plus
any it can borrow from threads left idle, which happens once the pool starts to
drain.  A task may borrow up to a share of the budget proportional to its
weight, such as its number of reads, relative to the heaviest task, so large
tasks get more of the idle threads.  Borrowed threads are returned when the
program finishes, so running tasks plus borrowed threads never exceed the
budget.

init() must be called in the parent before the pool is forked, so the workers
share the counters.
"""
import multiprocessing as mp
from contextlib import contextmanager

_threads = 1
_max_weight = 1
_state = None       # shared [running tasks, borrowed threads]
_task_weight = 1    # weight of the task running in this process


def init(threads, max_weight=1):
    "set the thread budget, and the weight of the heaviest task, for the next pool"
    global _threads, _max_weight, _state
    _threads = max(1, threads)
    _max_weight = max(1, max_weight)
    _state = mp.Array('i', [0, 0])


def task_start(weight=1):
    "record that a task of the given weight has started in this process"
    global _task_weight
    _task_weight = weight
    if _state is not None:
        with _state.get_lock():
            _state[0] += 1


def task_done():
    "record that the task running in this process has finished"
    if _state is not None:
        with _state.get_lock():
            _state[0] -= 1


@contextmanager
def program_threads():
    """Context manager giving the number of threads a program run by the current
    task may use: one, plus threads borrowed from the idle part of the budget, up
    to the task's weighted share.  The borrowed threads are returned on exit."""
    extra = 0
    if _state is not None:
        want = max(1, round(_threads * _task_weight / _max_weight))
        with _state.get_lock():
            extra = max(0, min(want - 1, _threads - _state[0] - _state[1]))
            _state[1] += extra
    try:
        yield 1 + extra
    finally:
        if extra:
            with _state.get_lock():
                _state[1] -= extra
