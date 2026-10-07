# Copyright 2006-2026 Mark Diekhans
"""Container index by sequence id, range, and optionally strand that
efficiently searches for overlapping entries.  Also contains code for
generating SQL where clauses to restrict by bin."""

##
# This code is a python reimplementation of binRange.{h,c} and hdb.c written
# by Jim Kent.
##

# from C code:
# There's a bin for each 128k (1<<17) segment. The next coarsest is 8x as big
# (1<<13).  That is for each 1M segment, for each 8M segment, for each 64M
# segment, for each 512M segment, and one top level bin for 4Gb.  Note, since
# start and end are int's, the practical limit is up to 2Gb-1, and thus, only
# four result bins on the second level. A range goes into the smallest bin it
# will fit in.

import sys
from collections import namedtuple
from flair.pycbio import PycbioException


class RangeFinderException(PycbioException):
    """Exception from RangeFinder"""
    pass

class _EntryNotFound(Exception):
    "low-level indication of entry not found"
    pass

class RemoveValueError(ValueError):
    "error when removed not found"
    def __init__(self, seqId, strand, entry):
        super(RemoveValueError, self).__init__(f"entry to remove not found: {seqId}:{entry.start}-{entry.end} strand={strand}, value: {repr(entry.value)}")


##
# Functions to translate ranges to bin numbers. Bin values are compatible with
# ones in UCSC database.
##

BIN_OFFSETS_BASIC = (512 + 64 + 8 + 1,
                     64 + 8 + 1,
                     8 + 1,
                     1, 0)

BIN_OFFSETS_EXTENDED = (4096 + 512 + 64 + 8 + 1,
                        512 + 64 + 8 + 1,
                        64 + 8 + 1, 8 + 1,
                        1, 0)

BIN_FIRST_SHIFT = 17  # How much to shift to get to finest bin.
BIN_NEXT_SHIFT = 3    # How much to shift to get to next larger bin.

BIN_BASIC_MAX_END = 512 * 1024 * 1024
BIN_OFFSET_TO_EXTENDED = 4681

def _calcBinForOffsets(start, end, baseOffset, offsets):
    "get the bin for a range"
    startBin = start >> BIN_FIRST_SHIFT
    endBin = (end - 1) >> BIN_FIRST_SHIFT
    for binOff in offsets:
        if (startBin == endBin):
            return baseOffset + binOff + startBin
        startBin >>= BIN_NEXT_SHIFT
        endBin >>= BIN_NEXT_SHIFT
    raise Exception("can't compute bin: start {}, end {} out of range".format(start, end))

def calcBin(start, end):
    "get the bin for a range"
    if end <= BIN_BASIC_MAX_END:
        return _calcBinForOffsets(start, end, 0, BIN_OFFSETS_BASIC)
    else:
        return _calcBinForOffsets(start, end, BIN_OFFSET_TO_EXTENDED, BIN_OFFSETS_EXTENDED)

def _getOverlappingBinsForOffsets(start, end, baseOffset, offsets):
    "generate bins for a range given a list of offsets"
    startBin = start >> BIN_FIRST_SHIFT
    endBin = (end - 1) >> BIN_FIRST_SHIFT
    for offset in offsets:
        yield (startBin + baseOffset + offset, endBin + baseOffset + offset)
        startBin >>= BIN_NEXT_SHIFT
        endBin >>= BIN_NEXT_SHIFT

def getOverlappingBins(start, end):
    """Generate bins for the range.  Each value is closed range of (startBin, endBin)"""
    if end <= BIN_BASIC_MAX_END:
        # contained in basic range
        for bins in _getOverlappingBinsForOffsets(start, end, 0, BIN_OFFSETS_BASIC):
            yield bins
        yield (BIN_OFFSET_TO_EXTENDED, BIN_OFFSET_TO_EXTENDED)
    else:
        if start < BIN_BASIC_MAX_END:
            # overlapping both basic and extended
            for bins in _getOverlappingBinsForOffsets(start, BIN_BASIC_MAX_END, 0, BIN_OFFSETS_BASIC):
                yield bins
        for bins in _getOverlappingBinsForOffsets(start, end, BIN_OFFSET_TO_EXTENDED, BIN_OFFSETS_EXTENDED):
            yield bins

def getOverlappingSqlExpr(binCol, seqCol, startCol, endCol, seq, start, end, *, placeHolder="?"):
    """generate an SQL expression for overlaps with the specified range, returning
    (expr, params).  The params must be bound to the query in the order returned;
    seq, start and end come from files and command lines, so they are never
    interpolated into the SQL.  In order to efficiently query, there should be an
    index on (seqName, bin).  Use placeHolder="%s" for a DB-API driver using
    format paramstyle, such as mysqlclient.
    """
    # bin numbers are computed here, not caller data, so they are literals
    parts = []
    for bins in getOverlappingBins(start, end):
        if bins[0] == bins[1]:
            parts.append("({}={})".format(binCol, bins[0]))
        else:
            parts.append("({}>={} and {}<={})".format(binCol, bins[0], binCol, bins[1]))
    expr = "(({}={}) and ({}<{}) and ({}>{}) and ({}))".format(seqCol, placeHolder,
                                                               startCol, placeHolder,
                                                               endCol, placeHolder,
                                                               " or ".join(parts))
    return expr, (seq, end, start)


class Entry(namedtuple("Entry",
                       ("start", "end", "value"))):
    "entry associating a range with a value"
    __slots__ = ()

    def overlaps(self, start, end):
        "test if the range is overlapped by the entry"
        return (end > self.start) and (start < self.end)

    def __str__(self):
        return str(self.start) + "-" + str(self.end) + ": " + str(self.value)


class RangeBins:
    """Range indexed container for a single sequence.  This using a binning
    scheme that implements spacial indexing. Based on UCSC hg browser binRange
    C module."""
    __slots__ = ("seqId", "strand", "minStart", "maxEnd", "buckets")

    def __init__(self, seqId, strand):
        self.seqId = seqId
        self.strand = strand
        self.minStart = sys.maxsize
        self.maxEnd = 0
        self.buckets = {}  # indexed by bin, each bucket is a list of Entry objects

    def add(self, entry):
        bin = calcBin(entry.start, entry.end)
        entries = self.buckets.get(bin)
        if entries is None:
            self.buckets[bin] = entries = []
        entries.append(entry)
        self.minStart = min(self.minStart, entry.start)
        self.maxEnd = max(self.maxEnd, entry.end)

    def _mergeEntries(self, newEntry, overEntries, mergeFunc):
        for e in overEntries:
            self.remove(e)
        entries = [newEntry] + overEntries
        return Entry(min((e.start for e in entries)),
                     max((e.end for e in entries)),
                     mergeFunc([e.value for e in entries]))

    def addMerge(self, entry, mergeFunc):
        overEntries = list(self.overlapping(entry.start, entry.end))  # must make copy
        if len(overEntries) > 0:
            entry = self._mergeEntries(entry, overEntries, mergeFunc)
        self.add(entry)

    def overlapping(self, start, end):
        "generator of entries overlapping the specified range"
        if (start < end):
            for bins in getOverlappingBins(start, end):
                for j in range(bins[0], bins[1] + 1):
                    bucket = self.buckets.get(j)
                    if bucket is not None:
                        for entry in bucket:
                            if entry.overlaps(start, end):
                                yield entry

    def remove(self, entry):
        """Remove an entry with the particular range and value"""
        try:
            bucket = self.buckets[(calcBin(entry.start, entry.end))]  # KeyError if no bucket
            bucket.remove(entry)  # ValueError if no value
        except (KeyError, IndexError, ValueError):
            raise _EntryNotFound()

    def values(self):
        "generator over all values"
        for bin in self.buckets.values():
            for entry in bin:
                yield entry.value

    def entries(self):
        "generator over all Entry objects"
        for bin in self.buckets.values():
            for entry in bin:
                yield entry

    def dump(self, fh):
        "print contents for debugging purposes"
        for bin in list(self.buckets.keys()):
            fh.write(self.seqId + " (" + str(self.strand) + ") bin=" + str(bin) + "\n")
            for entry in self.buckets[bin]:
                fh.write("\t" + str(entry) + "\n")


class RangeFinder:
    """Container index by sequence id, range, and optionally strand.
    All entries added to the object must either have strand or not
    have strand.  A query without strand will find all overlapping
    entries on either strand if strand was specified when adding entries.
    """
    validStrands = set((None, "+", "-"))

    def __init__(self, *, mergeFunc=None):
        "mergeFunc must take a list of values and returned a new value"
        self.mergeFunc = mergeFunc
        self.haveStrand = None
        self.seqBins = {}

    def _checkStrand(self, strand):
        if strand not in (None, "+", "-"):
            raise RangeFinderException(f"invalid strand: '{strand}'")

    @staticmethod
    def _binKey(seqId, strand):
        return (seqId, strand)

    def _getBin(self, seqId, strand):
        "return None if no bin for seq+strand"
        return self.seqBins.get(self._binKey(seqId, strand))

    def _getBinsForAdd(self, seqId, start, end, strand):
        "setup for adding entries"
        self._checkStrand(strand)
        if self.haveStrand is None:
            self.haveStrand = (strand is not None)
        elif self.haveStrand != (strand is not None):
            raise RangeFinderException("all RangeFinder entries must all either have strand or not have strand")
        key = self._binKey(seqId, strand)
        bins = self.seqBins.get(key)
        if bins is None:
            self.seqBins[key] = bins = RangeBins(seqId, strand)
        return bins

    def add(self, seqId, start, end, value, strand=None):
        "add an entry for a sequence and range, and optional strand"
        bins = self._getBinsForAdd(seqId, start, end, strand)
        bins.add(Entry(start, end, value))

    def addMerge(self, seqId, start, end, value, strand=None):
        """add an entry for a sequence and range, and optional strand,
        merging with existing overlapping entries"""
        if self.mergeFunc is None:
            raise RangeFinderException("RangeFinder does not have a mergeFunc, can't merge entries")
        bins = self._getBinsForAdd(seqId, start, end, strand)
        bins.addMerge(Entry(start, end, value), self.mergeFunc)

    def addCoords(self, coords, value, strand=None):
        """Added using Coords object.  Note that coords.strands is direction
        for sequence, not orientation of annotation on a strand"""
        if coords.strand == '-':
            coords = coords.reverse()
        self.add(coords.name, coords.start, coords.end, value, strand)

    def addMergeCoords(self, coords, value, strand=None):
        """add an entry using coords, and optional strand, merging with
        existing overlapping entries.  Note that coords.strands is direction for
        sequence, not orientation of annotation on a strand
        """
        if coords.strand == '-':
            coords = coords.reverse()
        self.addMerge(coords.name, coords.start, coords.end, value, strand)

    def _overlappingEntriesStrand(self, seqId, start, end, strand):
        "check overlap on specific strand, which might be None"
        bins = self.seqBins.get((seqId, strand))
        if bins is not None:
            for entry in bins.overlapping(start, end):
                yield entry

    def _overlappingEntriesBothStrands(self, seqId, start, end):
        "return range overlaps, checking both strands"
        yield from self._overlappingEntriesStrand(seqId, start, end, '+')
        yield from self._overlappingEntriesStrand(seqId, start, end, '-')

    def overlappingEntries(self, seqId, start, end, strand=None):
        """Return generator over Entry objects overlapping the specified range on
        seqId, optional strand"""
        self._checkStrand(strand)
        if self.haveStrand and (strand is None):
            # must check on both strands
            return self._overlappingEntriesBothStrands(seqId, start, end)
        else:
            # must only check a specifc strand, or no strand to check
            if not self.haveStrand:
                strand = None  # no strand to check
            return self._overlappingEntriesStrand(seqId, start, end, strand)

    def overlappingEntriesByCoords(self, coords, strand=None):
        """generator over Entry object overlapping using a Coords object"""
        return self.overlappingEntries(coords.name, coords.start, coords.end, strand)

    def overlapping(self, seqId, start, end, strand=None):
        """Return generator over values overlapping the specified range on
        seqId, optional strand"""
        for entry in self.overlappingEntries(seqId, start, end, strand=strand):
            yield entry.value

    def overlappingByCoords(self, coords, strand=None):
        """generator over values overlapping using a Coords object"""
        return self.overlapping(coords.name, coords.start, coords.end, strand=strand)

    def _removeFromSeqBin(self, seqId, strand, entry):
        # strand maybe None
        bins = self.seqBins.get((seqId, strand))
        if bins is None:
            raise _EntryNotFound()
        bins.remove(entry)

    def _removeBothStrands(self, seqId, entry):
        "remove an entry, checking both strands"
        try:
            self._removeFromSeqBin(seqId, '+', entry)
        except _EntryNotFound:
            self._removeFromSeqBin(seqId, '-', entry)

    def _remove(self, seqId, strand, entry):
        if self.haveStrand and (strand is None):
            # must check on both strands
            self._removeBothStrands(seqId, entry)
        else:
            # must only check a specifc strand, or no strand to check
            if not self.haveStrand:
                strand = None  # no strand to check
            self._removeFromSeqBin(seqId, strand, entry)

    def remove(self, seqId, start, end, value, strand=None):
        """remove an entry with the particular range and value, value error if not found"""
        self._checkStrand(strand)
        entry = Entry(start, end, value)
        try:
            self._remove(seqId, strand, entry)
        except _EntryNotFound:
            raise RemoveValueError(seqId, strand, entry)

    def values(self):
        "generator over all Entries object"
        for bins in self.seqBins.values():
            yield from bins.values()

    def entries(self):
        "generator over all Entriy objects"
        for bins in self.seqBins.values():
            yield from bins.entries()

    def getSeqIds(self):
        """get set of sequences"""
        return frozenset([k[0] for k in self.seqBins.keys()])

    def getSeqRange(self, seqId, strand=None):
        """Return the minimum start and maximum end for a given sequence.  Useful for
        look for closest entry by checking range.  If strand is not None and
        we have strands in bins, return just for that strand, otherwise for
        sequence regardless of strand.  Returns (None, None) if sequence is not
        present"""

        def _getRange(seqId, keyStrand):
            bins = self._getBin(seqId, keyStrand)
            if bins is None:
                return None, None
            else:
                return bins.minStart, bins.maxEnd

        self._checkStrand(strand)
        if not self.haveStrand:
            return _getRange(seqId, None)
        elif strand is None:
            # a sequence may have entries on only one strand
            r1 = _getRange(seqId, '+')
            r2 = _getRange(seqId, '-')
            if r1[0] is None:
                return r2
            elif r2[0] is None:
                return r1
            else:
                return (min(r1[0], r2[0]), max(r1[1], r2[1]))
        else:
            return _getRange(seqId, strand)

    def dump(self, fh):
        "print contents for debugging purposes"
        for bins in list(self.seqBins.values()):
            bins.dump(fh)
