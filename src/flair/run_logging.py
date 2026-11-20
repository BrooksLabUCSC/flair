"""What a run reports about itself: the format of a log record, what is logged at
the start of a run, and the level given to the flair programs a run starts.

The default level is INFO, so that a run which takes hours states what it is doing
without being asked.  A record is stamped with the time, and a warning or worse with
its level, so that a warning stands out in the progress of a long run.
"""

import logging
import os
import shlex
import sys
import time
import pipettor
from flair import VERSION
from flair.process_pool import fork_context

_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'


class FlairLogFormatter(logging.Formatter):
    """The time on every record, and the level on a warning, so that a warning stands
    out in the progress of a long run.  An INFO record is that progress, and the
    CRITICAL record that ends a run names its own error, so neither is marked."""

    def __init__(self):
        super().__init__(fmt='%(asctime)s %(message)s', datefmt=_TIME_FORMAT)
        self._marked = logging.Formatter(fmt='%(asctime)s %(levelname)s: %(message)s',
                                         datefmt=_TIME_FORMAT)

    def format(self, record):
        if logging.WARNING <= record.levelno < logging.CRITICAL:
            return self._marked.format(record)
        return super().format(record)


def debug_logging():
    "whether debug logging was asked for"
    return logging.getLogger().getEffectiveLevel() <= logging.DEBUG


def quiet_pipettor_failures():
    """pipettor logs a failed process at ERROR, with a traceback of the flair code
    that started it; run_program reports the failure as a tool error instead, so the
    record is dropped unless debug logging was asked for."""
    if not debug_logging():
        logging.getLogger(pipettor.LOGGER_NAME).setLevel(logging.CRITICAL)


def configure_logging():
    """Finish configuring logging, once the command line options have been applied:
    the flair format on each handler, and pipettor's own record of a failed process
    dropped.  Called by each program after its command line is parsed."""
    for handler in logging.getLogger().handlers:
        handler.setFormatter(FlairLogFormatter())
    quiet_pipettor_failures()


def log_run_start():
    """The version and the command line, which a long run's log, or a user reporting
    a problem with one, would otherwise state nothing about.  The command line is
    quoted as a shell would take it, so it can be run again as it is logged."""
    argv = [os.path.basename(sys.argv[0])] + sys.argv[1:]
    logging.info('FLAIR %s: %s', VERSION, shlex.join(argv))


def child_log_level_opt():
    """The --log-level option for a flair program this run starts, so that its
    progress and its diagnostics are those the run was asked for.  The level is
    passed as a number, which covers a level that has no name."""
    return f'--log-level={logging.getLogger().getEffectiveLevel()}'


def plural_noun(count, noun, other=None):
    "a noun in the number a count takes, other being its plural when adding s is wrong"
    return noun if count == 1 else (other or noun + 's')


def plural(count, noun, other=None):
    "a count and its noun, in the number the count takes"
    return f'{count} {plural_noun(count, noun, other)}'


# the tasks of the parallel stage in progress, counted in shared memory because each
# is reported by the worker that finishes it, and the workers are forked
_tasks_done = fork_context().Value('i', 0)
_tasks_total = fork_context().Value('i', 0)


# the reads a run has read, counted by the workers that read them
_reads_read = fork_context().Value('q', 0)


def count_reads_read(count):
    "add to the reads this run has read, from a worker or from the parent"
    if count:
        with _reads_read.get_lock():
            _reads_read.value += count


def report_reads_used(used, what):
    """How many of the reads a run read it reported, which is otherwise only in the
    outputs: a run that reports a quarter of its reads looks like one that reports
    all of them, and a bad library, a bad threshold and a bad annotation look
    alike."""
    read = _reads_read.value
    if read > 0:
        logging.info('of %s read, %d %s (%.0f%%)', plural(read, 'read'), used, what, 100 * used / read)


def start_task_progress(total):
    "the tasks of a parallel stage about to start, counted from one again"
    with _tasks_done.get_lock():
        _tasks_done.value = 0
        _tasks_total.value = total


def report_task_done(what, detail, start_time):
    """One line as a task of a parallel stage finishes, which for a run of hours is
    the only statement of how far it has got.  Called in the worker that ran the task,
    so the count of those finished is shared memory."""
    with _tasks_done.get_lock():
        _tasks_done.value += 1
        done = _tasks_done.value
    logging.info('%s %d of %d: %s, %.1f s', what, done, _tasks_total.value, detail,
                 time.time() - start_time)
