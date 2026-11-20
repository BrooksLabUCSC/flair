"""Stopping a run on a signal.

A run removes its temporary directory in a finally, which an error reaches and a
signal does not: a run stopped by a job scheduler, or at the terminal, left the
reads and alignments of every partition behind, under a name that states only which
output they belong to.

SIGTERM and SIGINT are raised as an error instead, so that the same finally runs.
SIGKILL cannot be taken, and a scheduler sends SIGTERM before it.
"""

import signal
from flair import FlairStopped


def _stop_run(signum, frame):
    raise FlairStopped(f'stopped by {signal.Signals(signum).name}')


def install_stop_handlers():
    "called by each program once its command line is parsed"
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, _stop_run)
