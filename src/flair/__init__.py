import os
from collections import namedtuple
from types import NoneType
from flair.pycbio import NoStackError

VERSION = "3.0.0+dev"

# Fixed minimum and maximums size for an intron sizes
# This could be configurable at some point, especially
# maximum size, do to crazy bad alignments
#
# Information glean from multiple sources,
#
# Minimal intron size
#   - smallest annotated human intron is 30bp in MST1L
#   - most indicate minimal size for most vertebrates
#     is longer than this.
#
# Maximum intron size
#   - below 10kb is common
#   - DOI:10.1371/journal.pone.0233978 finds some brain
#     expression of intron > 1mb
#
MIN_INTRON_SIZE = 2
MAX_INTRON_SIZE = 1000000

class FlairError(Exception):
    """General error condition in FLAIR"""
    pass

class FlairInputDataError(FlairError, NoStackError):
    """Error in FLAIR input data"""
    pass

class FlairNotImplementedError(FlairError, NoStackError):
    """An option is accepted by the argument parser but does nothing.  Raising is
    better than ignoring it silently, and better than dropping the option, which would
    make an existing command line fail with an argparse message that does not say the
    feature is missing."""
    pass

def set_unix_path():
    "add programs in package to PATH."

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explicit
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(os.path.realpath(__file__)) + ':' + os.environ["PATH"]


class PosRange(namedtuple("PosRange",
                          ("start", "end"))):
    """0-base, 1/2 open range of sequence positions"""

    def __new__(cls, start, end):
        assert start <= end  # allows zero length
        return super(PosRange, cls).__new__(cls, start, end)

    def __len__(self):
        return self.end - self.start


class SeqRange(namedtuple("SeqRange",
                          ("name", "start", "end", "strand"))):
    """Range within a sequence, including optional strand. Duck-type inheritance
    from PosRange"""

    def __new__(cls, name, start, end, strand=None):
        assert start <= end  # allows zero length
        assert isinstance(strand, (str, NoneType)) and ((strand is None) or (strand in ('+', '-', '.', None)))
        return super(SeqRange, cls).__new__(cls, name, start, end, strand)


def range_overlap(a, b):
    """Return the overlap length between two ranges with .start/.end attributes. Returns 0 if no overlap."""
    return max(0, min(a.end, b.end) - max(a.start, b.start))
