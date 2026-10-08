import multiprocessing as mp
from flair import thread_share


def test_busy_pool_gives_one_thread():
    thread_share.init(8, tasks=20)
    thread_share.task_start()
    with thread_share.program_threads() as n:
        assert n == 1


def test_idle_threads_are_borrowed_and_returned():
    thread_share.init(8, tasks=3)         # 5 threads idle
    thread_share.task_start()
    with thread_share.program_threads() as n:
        assert n == 6                     # its own thread plus the 5 idle ones
        with thread_share.program_threads() as m:
            assert m == 1                 # nothing left to borrow
    with thread_share.program_threads() as n:
        assert n == 6                     # returned on exit


def test_borrowing_is_proportional_to_weight():
    thread_share.init(8, tasks=2, max_weight=100)
    thread_share.task_start(25)           # a quarter of the heaviest task
    with thread_share.program_threads() as n:
        assert n == 2
    thread_share.task_start(100)          # the heaviest task can take them all
    with thread_share.program_threads() as n:
        assert n == 7                     # 8 minus the 2 unfinished tasks, plus its own


def test_finished_tasks_free_threads():
    thread_share.init(4, tasks=10)
    for _ in range(6):
        thread_share.task_done()
    with thread_share.program_threads() as n:
        assert n == 1                     # 4 tasks left keep all 4 workers busy
    thread_share.task_done()
    with thread_share.program_threads() as n:
        assert n == 2                     # 3 left: one worker idle


def test_single_thread():
    thread_share.init(1, tasks=1)
    thread_share.task_start()
    with thread_share.program_threads() as n:
        assert n == 1


# the startup race from the PR review: in a forked pool of 4 with 4 tasks, the
# first worker reaches its program before any other worker has started a task;
# it must not borrow the threads those workers are about to use
def _first_task():
    thread_share.task_start()
    try:
        with thread_share.program_threads() as n:
            return n
    finally:
        thread_share.task_done()


def test_no_borrowing_before_other_workers_start():
    thread_share.init(4, tasks=4)
    with mp.get_context('fork').Pool(4) as pool:
        assert pool.apply_async(_first_task).get(timeout=60) == 1
