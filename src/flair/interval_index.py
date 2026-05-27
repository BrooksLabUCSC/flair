"""
Lightweight per-chrom interval index backed by ruranges.

Supports the build-once/query-many pattern used in flair.  Intervals are
half-open ``[start, end)``; queries return the attached payload objects.
"""
import numpy as np
from ruranges import numpy as _rn

_COORD_DTYPE = np.int64


class IntervalIndex:
    __slots__ = ("_starts", "_ends", "_data",
                 "_arr_starts", "_arr_ends", "_dirty")

    def __init__(self):
        self._starts = []
        self._ends = []
        self._data = []
        self._arr_starts = None
        self._arr_ends = None
        self._dirty = False

    def __len__(self):
        return len(self._data)

    def add(self, start, end, data):
        self._starts.append(start)
        self._ends.append(end)
        self._data.append(data)
        self._dirty = True

    def _freeze(self):
        if not self._dirty and self._arr_starts is not None:
            return
        self._arr_starts = np.asarray(self._starts, dtype=_COORD_DTYPE)
        self._arr_ends = np.asarray(self._ends, dtype=_COORD_DTYPE)
        self._dirty = False

    def overlap(self, start, end, slack=0):
        """Return list of payloads whose interval overlaps ``[start, end)``.
        ``slack`` extends both sides of the query range."""
        if not self._data:
            return []
        self._freeze()
        qs = np.array([start], dtype=_COORD_DTYPE)
        qe = np.array([end], dtype=_COORD_DTYPE)
        idx1, _ = _rn.overlaps(
            starts=self._arr_starts, ends=self._arr_ends,
            starts2=qs, ends2=qe,
            slack=slack, sort_output=False)
        data = self._data
        return [data[i] for i in idx1]

    def items(self):
        "yield (start, end, data) for every interval"
        for start, end, data in zip(self._starts, self._ends, self._data):
            yield start, end, data
