from flair import thread_share


def _start_tasks(n, weight=1):
    for _ in range(n):
        thread_share.task_start(weight)


def test_busy_pool_gives_one_thread():
    thread_share.init(8)
    _start_tasks(8)
    with thread_share.program_threads() as n:
        assert n == 1


def test_idle_threads_are_borrowed_and_returned():
    thread_share.init(8)
    _start_tasks(3)                       # 5 threads idle
    with thread_share.program_threads() as n:
        assert n == 6                     # its own thread plus the 5 idle ones
        with thread_share.program_threads() as m:
            assert m == 1                 # nothing left to borrow
    with thread_share.program_threads() as n:
        assert n == 6                     # returned on exit


def test_borrowing_is_proportional_to_weight():
    thread_share.init(8, max_weight=100)
    _start_tasks(1, weight=25)            # a quarter of the heaviest task
    with thread_share.program_threads() as n:
        assert n == 2
    thread_share.task_start(100)          # the heaviest task can take them all
    with thread_share.program_threads() as n:
        assert n == 7                     # 8 minus the 2 running tasks, plus its own


def test_finished_tasks_free_threads():
    thread_share.init(4)
    _start_tasks(4)
    for _ in range(3):
        thread_share.task_done()
    with thread_share.program_threads() as n:
        assert n == 4


def test_single_thread():
    thread_share.init(1)
    _start_tasks(1)
    with thread_share.program_threads() as n:
        assert n == 1
