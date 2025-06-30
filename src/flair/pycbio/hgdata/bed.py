# Copyright 2006-2025 Mark Diekhans
import copy
from collections import deque, defaultdict, namedtuple
from flair.pycbio import PycbioException
from flair.pycbio.sys.color import Color
from flair.pycbio.tsv.tabFile import TabFile, TabFileReader
from flair.pycbio.hgdata.autoSql import intArraySplit, intArrayJoin

# FIXME: not complete, needs tests
# FIXME: really need a better way to deal with derived classes than extraCols


bed12Columns = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
                "thickEnd", "reserved", "blockCount", "blockSizes", "chromStarts")

class BedException(PycbioException):
    """Error parsing or operating on a BED"""
    pass

def defaultIfNone(v, dflt=""):
    "also converts to a string"
    return str(v) if v is not None else str(dflt)

def encodeRow(row):
    """convert a list of values to a list of strings, making None empty.
    """
    return [str(v) if v is not None else "" for v in row]

def _fmtItemRgb(itemRgb):
    "allows itemRgb to be a Color, None, a number, or a string"
    if itemRgb is None:
        return "0"
    elif isinstance(itemRgb, Color):
        return itemRgb.toRgb8Str()
    else:
        return itemRgb

class BedBlock(namedtuple("Block", ("start", "end"))):
    """A block in the BED.  Coordinates are absolute, not relative and are transformed on write.
    """
    __slots__ = ()

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return "{}-{}".format(self.start, self.end)

class Bed:
    """Object wrapper for a BED record.  ExtraCols is a vector of extra
    columns to add.  Special columns be added by extending and overriding
    parse() and toRow(), numColumns to do special handling.

    Columns maybe sparsely specified, with ones up to numStdCols defaulted.

    For BEDs with extra columns not handled by derived are stored in extraCols.
    If extra columns is a tuple, include namedtuple, it is stored as-is, otherwise
    a copy is stored.

    itemRgb can be a string or Color object
    """
    __slots__ = ("chrom", "chromStart", "chromEnd", "name", "score",
                 "strand", "thickStart", "thickEnd", "itemRgb", "blocks",
                 "extraCols", "numStdCols")

    def __init__(self, chrom, chromStart, chromEnd, name=None, *, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blocks=None, extraCols=None,
                 numStdCols=None):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blocks = copy.copy(blocks)
        self.extraCols = extraCols if isinstance(extraCols, tuple) else copy.copy(extraCols)
        self.numStdCols = self._calcNumStdCols(numStdCols)

    def _calcNumStdCols(self, specNumStdCols):
        # computer based on maximum specified
        if self.blocks is not None:
            numStdCols = 12
        elif self.itemRgb is not None:
            numStdCols = 9
        elif self.thickStart is not None:
            numStdCols = 8
        elif self.strand is not None:
            numStdCols = 6
        elif self.score is not None:
            numStdCols = 5
        elif self.name is not None:
            numStdCols = 4
        else:
            numStdCols = 3
        if specNumStdCols is not None:
            if not (3 <= specNumStdCols <= 12):
                raise BedException(f"numStdCols must be in the range 3 to 12, got {specNumStdCols}")
            if numStdCols > specNumStdCols:
                raise BedException(f"numStdCols was specified as {specNumStdCols}, however the arguments supplied require {numStdCols} standard columns")
            if specNumStdCols > numStdCols:
                numStdCols = specNumStdCols
        return numStdCols

    def addBlock(self, start, end):
        """add a new block"""
        assert start < end
        assert start >= self.chromStart
        assert end <= self.chromEnd
        blk = BedBlock(start, end)
        if self.blocks is None:
            self.blocks = []
            if self.numStdCols < 12:
                self.numStdCols = 12
        self.blocks.append(blk)
        return blk

    @property
    def numColumns(self):
        """Returns the number of columns in the BED when formatted as a row."""
        # exclude extraCols
        n = self.numStdCols
        if self.extraCols is not None:
            n += len(self.extraCols)
        return n

    def _getBlockColumns(self):
        relStarts = []
        sizes = []
        for blk in self.blocks:
            relStarts.append(str(blk.start - self.chromStart))
            sizes.append(str(len(blk)))
        return str(len(self.blocks)), intArrayJoin(sizes), intArrayJoin(relStarts)

    def _defaultBlockColumns(self):
        return "1", str(self.chromEnd - self.chromStart) + ',', "0,"

    def toRow(self):
        row = [self.chrom, str(self.chromStart), str(self.chromEnd)]
        if self.numStdCols >= 4:
            row.append(str(self.name) if self.name is not None else f"{self.chrom}:{self.chromStart}-{self.chromEnd}")
        if self.numStdCols >= 5:
            row.append(defaultIfNone(self.score, 0))
        if self.numStdCols >= 6:
            row.append(defaultIfNone(self.strand, '+'))
        if self.numStdCols >= 8:
            row.append(defaultIfNone(self.thickStart, self.chromEnd))
            row.append(defaultIfNone(self.thickEnd, self.chromEnd))
        if self.numStdCols >= 9:
            row.append(_fmtItemRgb(self.itemRgb))
        if self.numStdCols >= 10:
            row.extend(self._getBlockColumns() if self.blocks is not None else self._defaultBlockColumns())
        if self.extraCols is not None:
            row.extend(encodeRow(self.extraCols))
        return row

    @staticmethod
    def _parseBlockColumns(chromStart, row):
        sizes = intArraySplit(row[10])
        relStarts = intArraySplit(row[11])
        blocks = []
        for i in range(len(relStarts)):
            start = chromStart + relStarts[i]
            blocks.append(BedBlock(start, start + sizes[i]))
        return blocks

    @classmethod
    def _parse(cls, row, numStdCols=None, *, fixScores=False):
        assert (numStdCols is None) or (3 <= numStdCols <= 12)
        if numStdCols is None:
            numStdCols = min(len(row), 12)
        if len(row) < numStdCols:
            raise BedException("expected at least {} columns, found {}: ".format(numStdCols, len(row)))
        chrom = row[0]
        chromStart = int(row[1])
        chromEnd = int(row[2])
        if numStdCols > 3:
            name = row[3]
        else:
            name = None
        if numStdCols > 4:
            try:
                # match browser behavior of converting to ints
                score = int(float(row[4]))
            except ValueError:
                if fixScores:
                    score = 0
                else:
                    raise
        else:
            score = None
        if numStdCols > 5:
            strand = row[5]
        else:
            strand = None
        if numStdCols > 7:
            thickStart = int(row[6])
            thickEnd = int(row[7])
        else:
            thickStart = None
            thickEnd = None
        if numStdCols > 8:
            itemRgb = row[8]
        else:
            itemRgb = None
        if numStdCols > 11:
            blocks = Bed._parseBlockColumns(chromStart, row)
        else:
            blocks = None
        if len(row) > numStdCols:
            extraCols = row[numStdCols:]
        else:
            extraCols = None
        return cls(chrom, chromStart, chromEnd, name=name, score=score, strand=strand,
                   thickStart=thickStart, thickEnd=thickEnd, itemRgb=itemRgb, blocks=blocks,
                   extraCols=extraCols, numStdCols=numStdCols)

    @classmethod
    def parse(cls, row, numStdCols=None, *, fixScores=False):
        """Parse a list of BED columns, as strings, into a Bed object.  If
        self.numStdCols is specified, only those columns are parsed and the
        remainder goes into extraCols.  Floating point scores are converted to
        ints to match UCSC browser behavior. If fixScores is True, non-numeric
        scores are converted to zero rather than generating an error."""
        try:
            return cls._parse(row, numStdCols=numStdCols, fixScores=fixScores)
        except Exception as ex:
            raise BedException(f"parsing of BED row failed: {row}") from ex

    def __str__(self):
        "return BED as a tab-separated string"
        return "\t".join(self.toRow())

    @property
    def start(self):
        """alias for chromStart"""
        return self.chromStart

    @property
    def end(self):
        """alias for chromEnd"""
        return self.chromEnd

    @property
    def blockCount(self):
        return len(self.blocks)

    @property
    def span(self):
        "distance from start to end"
        return self.chromEnd - self.chromStart

    @property
    def coverage(self):
        """number of bases covered"""
        if self.blocks is None:
            return self.span
        else:
            return sum([len(b) for b in self.blocks])

    def getGaps(self):
        """return a tuple of BedBlocks for the coordinates of the gaps between
        the blocks, which are often introns"""
        gaps = []
        prevBlk = None
        for blk in self.blocks:
            if prevBlk is not None:
                gaps.append(BedBlock(prevBlk.end, blk.start))
            prevBlk = blk
        return tuple(gaps)

    def write(self, fh):
        """write BED to a tab-separated file"""
        fh.write(str(self))
        fh.write('\n')

    @staticmethod
    def genome_sort_key(bed):
        return bed.chrom, bed.chromStart

def BedReader(fspec, numStdCols=None, bedClass=Bed, *, fixScores=False):
    """Generator to read BED objects loaded from a tab-file or file-like
    object.  See Bed.parse()."""
    for bed in TabFileReader(fspec, rowClass=lambda r: bedClass.parse(r, numStdCols=numStdCols, fixScores=fixScores),
                             hashAreComments=True, skipBlankLines=True):
        yield bed


class BedTable(TabFile):
    """Table of BED objects loaded from a tab-file
    """

    def _mkNameIdx(self):
        self.nameMap = defaultdict(list)
        for bed in self:
            self.nameMap[bed.name].append(bed)

    def __init__(self, fileName, nameIdx=False, numStdCols=None):
        super(BedTable, self).__init__(fileName, rowClass=lambda r: Bed.parse(r, numStdCols=numStdCols),
                                       hashAreComments=True, skipBlankLines=True)
        self.nameMap = None
        if nameIdx:
            self._mkNameIdx()

    def getByName(self, name):
        "get *tuple* of BEDs by name"
        if name in self.nameMap:
            return tuple(self.nameMap[name])
        else:
            return ()

def bedFromPsl(psl, *, extraCols=None):
    "create a BED12 from PSL, optionally adding extra columns"
    if psl.tStrand == '-':
        psl = psl.reverseComplement()
    blks = [BedBlock(pb.tStart, pb.tEnd) for pb in psl.blocks]

    return Bed(psl.tName, psl.tStart, psl.tEnd, name=psl.qName, strand=psl.qStrand,
               blocks=blks, extraCols=extraCols, numStdCols=12)

###
# bedMergeBlocks
###
def _bedMergeCheckCompat(bed0, bed):
    if bed.numStdCols != bed0.numStdCols:
        raise BedException(f"attempt to merge BEDs number of standard columns: {bed.name}: {bed.numStdCols} != {bed0.name} {bed0.numStdCols}")
    if bed.chrom != bed0.chrom:
        raise BedException(f"attempt to merge BEDs on different chromosomes: {bed.name}: {bed.chrom} != {bed0.name} {bed0.chrom}")
    if bed.strand != bed0.strand:
        raise BedException(f"attempt to merge BEDs on different strands: {bed.name}: {bed.strand} != {bed0.name} {bed0.strand}")

def _bedMergeBuildCursors(beds):
    """Build a vector of cursors and validated sequence and strand compatibility"""
    cursors = []
    bed0 = beds[0]
    minStart = bed0.chromStart
    maxEnd = bed0.chromEnd
    for bed in beds:
        _bedMergeCheckCompat(bed0, bed)
        cursors.append(deque(bed.blocks))
        minStart = min(bed.chromStart, minStart)
        maxEnd = max(bed.chromEnd, maxEnd)
    return cursors, minStart, maxEnd

def _bedMergeFindNextStart(cursors):
    """locate the lowest starting position, or None if all cursors are empty"""
    nextStart = None
    for cursor in cursors:
        if (len(cursor) > 0) and ((nextStart is None) or (cursor[0].start < nextStart)):
            nextStart = cursor[0].start
    return nextStart

def _bedMergePass(cursors, mergeBlkEnd):
    """make on pass over cursors, seeing if range can be updated with
    overlapping or adjacent blocks.  Since processed blocks are removed, from
    cursor, the block at the head of the cursor must either overlap or be
    after the block.  Multiple passed handle transitive joins.
    """
    updated = False
    for cursor in cursors:
        if (len(cursor) > 0) and (cursor[0].start <= mergeBlkEnd):
            mergeBlkEnd = max(cursor[0].end, mergeBlkEnd)
            cursor.popleft()
            updated = True
    return mergeBlkEnd, updated

def _bedMergeBuildBlk(cursors, mergeBlkStart):
    "build one block"
    mergeBlkEnd = mergeBlkStart + 1
    while True:
        mergeBlkEnd, updated = _bedMergePass(cursors, mergeBlkEnd)
        if not updated:
            break
    return BedBlock(mergeBlkStart, mergeBlkEnd)

def _bedMergeBuildBlks(cursors):
    "build all blocks, run untils cursors are empty"
    mergedBlocks = []
    nextStart = _bedMergeFindNextStart(cursors)
    while nextStart is not None:
        mergedBlocks.append(_bedMergeBuildBlk(cursors, nextStart))
        nextStart = _bedMergeFindNextStart(cursors)
    return mergedBlocks

def _bedMergeThickBounds(beds, maxEnd):
    thickStart = thickEnd = None
    for bed in beds:
        if (bed.thickStart is not None) and (bed.thickStart < bed.thickEnd):
            if thickStart is None:
                thickStart, thickEnd = bed.thickStart, bed.thickEnd
            else:
                thickStart, thickEnd = min(thickStart, bed.thickStart), max(thickEnd, bed.thickEnd)
    if thickStart is None:
        thickStart = thickEnd = maxEnd
    return thickStart, thickEnd

def bedMergeBlocks(name, beds):
    """Merge blocks from multiple BED into a single BED.  The BEDs must be on
    the same sequence and strand.  The thickStart and thickStop will be the
    minimum and maximum seen, unless they are zero length, in which case, it
    will remain zero length.  Useful with RangeFinder merge functionality."""
    if len(beds) == 0:
        raise BedException("can't merge list with no BEDs")
    cursors, minStart, maxEnd = _bedMergeBuildCursors(beds)
    mergedBlocks = _bedMergeBuildBlks(cursors)
    thickStart, thickEnd = _bedMergeThickBounds(beds, maxEnd)

    bed0 = beds[0]
    return Bed(bed0.chrom, minStart, maxEnd, name=name, strand=bed0.strand,
               thickStart=thickStart, thickEnd=thickEnd, itemRgb=bed0.itemRgb,
               blocks=mergedBlocks, numStdCols=12)
