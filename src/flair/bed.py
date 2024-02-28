# based on bed.py by Mark Diekhans at https://github.com/diekhans/pycbio.git
import copy
from collections import defaultdict, namedtuple

bed12Columns = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
				"thickEnd", "reserved", "blockCount", "blockSizes", "chromStarts")

def intArraySplit(commaStr):
	"Splits a comma separated string into a list of ints"
	return([int(x) for x in commaStr.rstrip(',').split(',')])

def intArrayJoin(ints):
	"""Formats a list of integers into a comma-separated string."""
	if ints is not None:
		strs = [str(i) for i in ints]
		return ",".join(strs) + ","
	else:
		return ","

def defaultIfNone(v, dflt=""):
	# also converts to a string
	return str(v) if v is not None else str(dflt)

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
	For BEDs with extra columns not handled by derived are stored in extraCols
	"""
	__slots__ = ("chrom", "chromStart", "chromEnd", "name", "score",
				 "strand", "thickStart", "thickEnd", "itemRgb", "blocks",
				 "extraCols", "numStdCols")

	def __init__(self, chrom, chromStart, chromEnd, name=None, *, score=None, strand=None,
				 thickStart=None, thickEnd=None, itemRgb='0', blocks=None, extraCols=None,
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
		self.extraCols = copy.copy(extraCols)

		if numStdCols is None:
			# compute based on what is specified
			numStdCols = 3
			if name is not None:
				numStdCols += 1
			if score is not None:
				numStdCols += 1
			if strand is not None:
				numStdCols += 1
			if thickStart is not None:
				numStdCols += 2
			if thickStart is not None:
				numStdCols += 2
			if itemRgb is not None:
				numStdCols += 1
			if blocks is not None:
				numStdCols += 3
		self.numStdCols = numStdCols

	def addBlock(self, start, end):
		"""add a new block"""
		assert start < end
		assert start >= self.chromStart
		assert end <= self.chromEnd
		blk = BedBlock(start, end)
		if self.blocks is None:
			self.blocks = []
		self.blocks.append(blk)
		return blk

	@property
	def numColumns(self):
		"""Returns the number of columns in the BED when formatted as a row."""
		# exclude extraCols
		n = self.numStdCols
		if n >= 9:
			n += 2  # blocks take up three columns
		if self.extraCols is not None:
			n += len(self.extraCols)
		return n

	def _getBlockColumns(self):
		relStarts = []
		sizes = []
		for blk in self.blocks:
			relStarts.append(str(blk.start - self.chromStart))
			sizes.append(str(len(blk)))
		return [intArrayJoin(sizes), intArrayJoin(relStarts)]

	def toRow(self):
		row = [self.chrom, str(self.chromStart), str(self.chromEnd)]
		if self.numStdCols >= 4:
			row.append(defaultIfNone(self.name))
		if self.numStdCols >= 5:
			row.append(defaultIfNone(self.score, 0))
		if self.numStdCols >= 6:
			row.append(defaultIfNone(self.strand, '+'))
		if self.numStdCols >= 7:
			row.append(defaultIfNone(self.thickStart, self.chromEnd))
		if self.numStdCols >= 8:
			row.append(defaultIfNone(self.thickEnd, self.chromEnd))
		if self.numStdCols >= 9:
			row.append(str(self.itemRgb))
		if self.numStdCols >= 10:
			row.append(str(len(self.blocks)))
			row.extend(self._getBlockColumns())
		if self.extraCols is not None:
			row.extend([defaultIfNone(c) for c in self.extraCols])
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
	def _parse(cls, row, numStdCols=None):
		assert (numStdCols is None) or (3 <= numStdCols <= 12)
		if numStdCols is None:
			numStdCols = min(len(row), 12)
		if len(row) < numStdCols:
			raise Exception("expected at least {} columns, found {}: ".format(numStdCols, len(row)))
		chrom = row[0]
		chromStart = int(row[1])
		chromEnd = int(row[2])
		if numStdCols > 3:
			name = row[3]
		else:
			name = None
		if numStdCols > 4:
			score = int(row[4])
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
	def parse(cls, row, numStdCols=None):
		"""Parse bed string columns into a bed object.  If self.numStdCols
		is specified, only those columns are parse and the remained goes
		to extraCols."""
		try:
			return cls._parse(row, numStdCols=numStdCols)
		except Exception as ex:
			raise Exception(f"parsing of BED row failed: {row}") from ex

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

	def write(self, fh):
		"""write BED to a tab-separated file"""
		fh.write(str(self))
		fh.write('\n')

	@staticmethod
	def genome_sort_key(bed):
		return bed.chrom, bed.chromStart

#def BedReader(fspec, numStdCols=None, bedClass=Bed):
#	"""Generator to read BED objects loaded from a tab-file or file-like
#	object.  See Bed.parse()."""
#	for bed in TabFileReader(fspec, rowClass=lambda r: bedClass.parse(r, numStdCols=numStdCols),
#							 hashAreComments=True, skipBlankLines=True):
#		yield bed


#class BedTable(TabFile):
#	"""Table of BED objects loaded from a tab-file
#	"""
#
#	def _mkNameIdx(self):
#		self.nameMap = defaultdict(list)
#		for bed in self:
#			self.nameMap[bed.name].append(bed)
#
#	def __init__(self, fileName, nameIdx=False, numStdCols=None):
#		super(BedTable, self).__init__(fileName, rowClass=lambda r: Bed.parse(r, numStdCols=numStdCols),
#									   hashAreComments=True, skipBlankLines=True)
#		self.nameMap = None
#		if nameIdx:
#			self._mkNameIdx()
#
#	def getByName(self, name):
#		"get *tuple* of BEDs by name"
#		if name in self.nameMap:
#			return tuple(self.nameMap[name])
#		else:
#			return ()
