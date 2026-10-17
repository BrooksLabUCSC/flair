"""
Lightweight per-chrom interval index with simple searching.

Supports the build-once/query-many pattern used in flair.  Intervals are
half-open ``[start, end)``; queries return the attached payload objects.

Coordinates are kept in arrays, and searched through a start-sorted order with
bisection, so a query only reads the intervals near it.  Besides being faster
than a scan, this matters for forked workers sharing a large index built by the
parent: reading array buffers doesn't touch Python objects, so their memory
pages aren't copied into each worker.  Call build() in the parent before
forking, so the sorted order is shared too.
"""
from array import array
from bisect import bisect_left, bisect_right


class IntervalIndex:
    __slots__ = ("_starts", "_ends", "_data",
                 "_order", "_sorted_starts", "_max_len", "_dirty")

    def __init__(self):
        self._starts = array('q')
        self._ends = array('q')
        self._data = []
        self._order = None          # indexes of the intervals, sorted by start
        self._sorted_starts = None  # starts in that order
        self._max_len = 0
        self._dirty = False

    def __len__(self):
        return len(self._data)

    def add(self, start, end, data):
        self._starts.append(start)
        self._ends.append(end)
        self._data.append(data)
        self._dirty = True

    def build(self):
        "build the search order now, rather than on the next query"
        order = sorted(range(len(self._data)), key=self._starts.__getitem__)
        self._order = array('q', order)
        self._sorted_starts = array('q', (self._starts[i] for i in order))
        self._max_len = max((e - s for s, e in zip(self._starts, self._ends)), default=0)
        self._dirty = False

    def overlap(self, start, end, slack=0):
        """Return list of payloads whose interval overlaps ``[start, end)``,
        in the order they were added.  ``slack`` extends both sides of the query
        range."""
        if not self._data:
            return []
        if self._dirty or self._order is None:
            self.build()
        # slack on the query, not on the intersection: + slack on the min widened
        # whichever end happened to be smaller, so an interval to the left of the
        # query was never reached while one to the right was
        qstart, qend = start - slack, end + slack
        # an overlapping interval starts before the query ends, and, being at most
        # _max_len long, after qstart - _max_len
        lo = bisect_right(self._sorted_starts, qstart - self._max_len)
        hi = bisect_left(self._sorted_starts, qend)
        hits = []
        for i in range(lo, hi):
            j = self._order[i]
            if max(self._starts[j], qstart) < min(self._ends[j], qend):
                hits.append(j)
        hits.sort()
        return [self._data[j] for j in hits]

    def items(self):
        "yield (start, end, data) for every interval"
        for start, end, data in zip(self._starts, self._ends, self._data):
            yield start, end, data
