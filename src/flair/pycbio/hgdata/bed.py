# Copyright 2006-2025 Mark Diekhans
import copy
from collections import defaultdict, namedtuple
from flair.pycbio import PycbioException
from flair.pycbio.sys.color import Color
from flair.pycbio.tsv.tabFile import TabFile, TabFileReader
from flair.pycbio.hgdata.autoSql import intArraySplit, intArrayJoin, strArrayJoin

# FIXME: not complete, needs tests
# FIXME: really need a better way to deal with derived classes than extraCols
# FIXME: FLAIR uses is an example of issues


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
    erow = []
    for v in row:
        if v is None:
            erow.append("")
        elif isinstance(v, (list, tuple)):
            erow.append(strArrayJoin(v))
        else:
            erow.append(str(v))
    return erow

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
    If extra columns is a tuple, it is stored as-is, otherwise a copy is stored.
    Non-string fields are converted to string on output. None is converted to
    an empty string.

    itemRgb can be a string or Color object
    """
    __slots__ = ("chrom", "chromStart", "chromEnd", "name", "score",
                 "strand", "thickStart", "thickEnd", "itemRgb", "blocks",
                 "extraCols", "numStdCols")

    def __init__(self, chrom, chromStart, chromEnd, name=None, *, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blocks=None, extraCols=(),
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
        self.extraCols = tuple(extraCols)  # copies unless it is already a tuple
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
        return self.numStdCols + len(self.extraCols)

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
        if len(self.extraCols) > 0:
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
    def _parse(cls, row, numStdCols=None, *, fixScores=False, skipExtraCols=False):
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
        if (not skipExtraCols) and (len(row) > numStdCols):
            extraCols = row[numStdCols:]
        else:
            extraCols = ()
        return cls(chrom, chromStart, chromEnd, name=name, score=score, strand=strand,
                   thickStart=thickStart, thickEnd=thickEnd, itemRgb=itemRgb, blocks=blocks,
                   extraCols=extraCols, numStdCols=numStdCols)

    @classmethod
    def parse(cls, row, numStdCols=None, *, fixScores=False, skipExtraCols=False):
        """Parse a list of BED columns, as strings, into a Bed object.  If
        self.numStdCols is specified, only those columns are parsed and the
        remainder goes into extraCols.  Floating point scores are converted to
        ints to match UCSC browser behavior. If fixScores is True, non-numeric
        scores are converted to zero rather than generating an error."""
        try:
            return cls._parse(row, numStdCols=numStdCols, fixScores=fixScores, skipExtraCols=skipExtraCols)
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
        if self.blocks is None:
            return 0
        else:
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

    def addExtraCols(self, cols):
        """append extra column values"""
        self.extraCols += tuple(cols)

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

def bedFromPsl(psl, *, extraCols=()):
    "create a BED12 from PSL, optionally adding extra columns"
    if psl.tStrand == '-':
        psl = psl.reverseComplement()
    blks = [BedBlock(pb.tStart, pb.tEnd) for pb in psl.blocks]

    return Bed(psl.tName, psl.tStart, psl.tEnd, name=psl.qName, strand=psl.qStrand,
               blocks=blks, extraCols=extraCols, numStdCols=12)

###
# bedMergeBlocks
###
def _bedMergeCollectBlocks(beds, stranded):
    """Collect all blocks from BEDs and validate compatibility."""
    bed0 = beds[0]
    all_blocks = list(bed0.blocks)
    for bed in beds[1:]:
        if bed.numStdCols != bed0.numStdCols:
            raise BedException(f"attempt to merge BEDs number of standard columns: {bed.name}: {bed.numStdCols} != {bed0.name} {bed0.numStdCols}")
        if bed.chrom != bed0.chrom:
            raise BedException(f"attempt to merge BEDs on different chromosomes: {bed.name}: {bed.chrom} != {bed0.name} {bed0.chrom}")
        if stranded and (bed.strand != bed0.strand):
            raise BedException(f"attempt to merge BEDs on different strands: {bed.name}: {bed.strand} != {bed0.name} {bed0.strand}")
        all_blocks.extend(bed.blocks)
    return all_blocks

def _bedMergeBlocks(all_blocks):
    """Sort and merge overlapping/adjacent blocks"""
    all_blocks.sort(key=lambda b: b.start)
    merged = []
    for blk in all_blocks:
        if merged and blk.start <= merged[-1].end:
            merged[-1] = BedBlock(merged[-1].start, max(merged[-1].end, blk.end))
        else:
            merged.append(blk)
    return merged

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

def bedMergeBlocks(name, beds, *, stranded=True):
    """Merge blocks from multiple BED into a single BED.  The BEDs must be on
    the same sequence.  If stranded is True (default), BEDs must be on the same
    strand; if False, different strands are allowed and the result strand is '+'.
    The thickStart and thickStop will be the minimum and maximum seen, unless
    they are zero length, in which case, it will remain zero length.
    Useful with RangeFinder merge functionality."""
    if len(beds) == 0:
        raise BedException("can't merge list with no BEDs")
    all_blocks = _bedMergeCollectBlocks(beds, stranded)
    mergedBlocks = _bedMergeBlocks(all_blocks)
    chromStart, chromEnd = mergedBlocks[0].start, mergedBlocks[-1].end
    thickStart, thickEnd = _bedMergeThickBounds(beds, chromEnd)

    bed0 = beds[0]
    strand = bed0.strand if stranded else '+'
    return Bed(bed0.chrom, chromStart, chromEnd, name=name, strand=strand,
               thickStart=thickStart, thickEnd=thickEnd, itemRgb=bed0.itemRgb,
               blocks=mergedBlocks, numStdCols=12)
