"""
Lightweight per-chrom interval index with simple searching.

Supports the build-once/query-many pattern used in flair.  Intervals are
half-open ``[start, end)``; queries return the attached payload objects.
"""


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

    def overlap(self, start, end, slack=0):
        """Return list of payloads whose interval overlaps ``[start, end)``.
        ``slack`` extends both sides of the query range."""
        if not self._data:
            return []
        overlaps = []
        # slack on the query, not on the intersection: + slack on the min widened
        # whichever end happened to be smaller, so an interval to the left of the
        # query was never reached while one to the right was
        for i in range(len(self._data)):
            if max(self._starts[i], start - slack) < min(self._ends[i], end + slack):
                overlaps.append(self._data[i])
        return overlaps

    def items(self):
        "yield (start, end, data) for every interval"
        for start, end, data in zip(self._starts, self._ends, self._data):
            yield start, end, data
