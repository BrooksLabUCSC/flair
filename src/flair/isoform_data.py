"""Read/exon structures and junction utilities, shared across FLAIR modules."""

from collections import namedtuple
from flair import PosRange
from flair.pycbio.hgdata.bed import Bed
from statistics import median
import pysam


POLYA_MIN_FRAC = 0.6
POLYA_SEARCH_WINDOW = 5
POLYA_MIN_LEN = 10
INTPRIM_MIN_FRAC = 0.6
INTPRIM_MIN_AS = 8
INTPRIM_SEARCH_WINDOW = 50

####
# basic types
####
class Junc(PosRange):
    """Stores start, end, just adds a type name to SeqRange for clearer code and error messages"""
    pass

class Exon(namedtuple("Exon", ("start", "end", "name"))):
    def __new__(cls, start, end, name=None):
        assert start <= end
        if name is None:
            name = ''
        return super(Exon, cls).__new__(cls, start, end, name)

    def __len__(self):
        return self.end - self.start


ISO_SRC_ANNOT = 'annot'
ISO_SRC_NOVEL = 'novel'


# class IsoIdSrc(namedtuple("IsoIdSrc",
#                           ("id", "src"))):
#     """isoform identifier along with the source of the isoform"""
#     # FIXME: it is unclear if this is the best way to store the information,
#     # this was create as a transition from iso (id) or (iso_id) (marker, id)
#     pass


def exons_to_juncs(exons):
    """Convert exon ranges to junctions"""
    return [Junc(exons[i].end, exons[i + 1].start)
            for i in range(len(exons) - 1)]


def bed_to_junctions(bed):
    # FIXME: a junctions object might be good
    return [Junc(bed.blocks[i - 1].end, bed.blocks[i].start)
            for i in range(1, bed.blockCount)]


def get_rgb(strand, junclen):
    if junclen == 0:
        return "99,99,99"
    elif strand == '+':
        return "27,158,119"
    else:
        return "217,95,2"


def get_bed_exons_from_juncs(juncs, start, end):
    if len(juncs) == 0:
        exon_starts = [0]
        exon_sizes = [end - start]
    else:
        exon_starts = [0] + [j.end - start for j in juncs]
        exon_sizes = ([juncs[0].start - start] + [juncs[i + 1].start - juncs[i].end
                                                  for i in range(len(juncs) - 1)] +
                      [end - juncs[-1].end])
    return exon_starts, exon_sizes


def get_bed_exons_from_exons(exons, start):
    exon_starts = [e.end - start for e in exons]
    exon_sizes = [e.end - e.start for e in exons]
    return exon_starts, exon_sizes


def get_reverse_complement(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S',
                'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    seq = seq.upper()
    new_seq = []
    for base in seq:
        new_seq.append(compbase[base])
    return ''.join(new_seq[::-1])


def get_sequence_for_exons(genome, chrom, strand, exons):
    trans_seq = ''.join([genome.fetch(chrom, e.start, e.end)
                         for e in exons])
    # FIXME: this upper cases only if reverse strand
    if strand == '-':
        trans_seq = get_reverse_complement(trans_seq)
    return trans_seq


def binary_search(query, data):
    """ Query is a coordinate interval. Binary search for the query in sorted data,
        which is a list of coordinates. Finishes when an overlapping value of query and
        data exists and returns the index in data. """
    # FIXME: uses python bisect module
    i = int(round(len(data) / 2))  # binary search prep
    lower, upper = 0, len(data)
    while True:
        if upper - lower < 2:  # stop condition but not necessarily found
            break
        if data[i].end < query.start:
            lower = i
            i = int(round((i + upper) / 2))
        elif data[i].start > query.end:
            upper = i
            i = int(round((lower + i) / 2))
        else:  # found
            break
    return i


def convert_to_bed(readrec):
    """Create and return a Bed object."""
    exon_starts, exon_sizes = get_bed_exons_from_juncs(readrec.juncs, readrec.start, readrec.end)
    bed = Bed(readrec.chrom, readrec.start, readrec.end, readrec.name,
              score=readrec.score, strand=readrec.strand,
              thickStart=readrec.start, thickEnd=readrec.end,
              itemRgb=get_rgb(readrec.strand, len(readrec.juncs)))
    for i in range(len(exon_starts)):
        blk_start = readrec.start + exon_starts[i]
        bed.addBlock(blk_start, blk_start + exon_sizes[i])
    return bed

def get_exons(readrec):
    """Return exons as list of Exon objects, computed from start, end, and juncs."""
    if not readrec.juncs:
        return [Exon(readrec.start, readrec.end)]
    exons = [Exon(readrec.start, readrec.juncs[0].start)]
    for i in range(len(readrec.juncs) - 1):
        exons.append(Exon(readrec.juncs[i].end, readrec.juncs[i + 1].start))
    exons.append(Exon(readrec.juncs[-1].end, readrec.end))
    return exons


# class JuncChain:
#     # keep on one copy of junction chain
#     _juncs_cache = {}

#     @classmethod
#     def _intern_juncs(cls, juncs):
#         return cls._juncs_cache.setdefault(juncs, juncs)

#     def __init__(self, chrom, strand, juncs):
#         self.chrom = chrom
#         self.strand = strand
#         self.juncs = self._intern_juncs(juncs)

def check_intprim(end_seq):
    i = 10
    while i < len(end_seq) and end_seq[:i].count('A') / i >= INTPRIM_MIN_FRAC:
        i += 1
    j = end_seq[:i].count('A')
    if j < INTPRIM_MIN_AS or j / i < INTPRIM_MIN_FRAC:
        return 0
    else:
        return j

def check_polyA(end_seq):
    if len(end_seq) < POLYA_SEARCH_WINDOW:
        return 0
    else:
        i = POLYA_SEARCH_WINDOW
        while i < len(end_seq) and end_seq[i - POLYA_SEARCH_WINDOW:i].count('A') / POLYA_SEARCH_WINDOW >= POLYA_MIN_FRAC:  # check rolling average of 5bp
            i += 1
        # return i
        j = end_seq[:i].rfind('A') + 1
        if j < POLYA_MIN_LEN:
            return 0
        else:
            return j


class ReadRec:
    """Read alignment with location, junction, and metadata fields.

    Stores chrom, start, end, name, score, strand, and juncs directly.
    Exons are computed on the fly from start, end, and juncs.
    A Bed record can be produced on demand via to_bed().

    Junction tuples are interned via a class-level cache so identical junction
    chains share a single tuple object.
    """

    # keep on one copy of junction chain
    _juncs_cache = {}

    @classmethod
    def _intern_juncs(cls, juncs):
        return cls._juncs_cache.setdefault(juncs, juncs)

    def __init__(self, chrom, strand, juncs, start, end, name, *, score=None, polyA=None, intprim=None):
        self.chrom = chrom
        self.strand = strand
        self.juncs = self._intern_juncs(juncs)
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.polyA = polyA  # (left int, right int)
        self.intprim = intprim  # (left int, right int)

    @property
    def exons(self):
        return get_exons(self)

    def reset_from_exons(self, exons):
        """Update ReadRec from a list of Exon objects."""
        self.start = exons[0].start
        self.end = exons[-1].end
        self.juncs = tuple(exons_to_juncs(sorted(exons)))

    def _get_both_intprim(read, genome):
        left_intprim, right_intprim = 0, 0
        if read.reference_start > INTPRIM_SEARCH_WINDOW:
            end_seq = get_reverse_complement(genome.fetch(read.reference_name, read.reference_start - INTPRIM_SEARCH_WINDOW, read.reference_start))
            left_intprim = check_intprim(end_seq)
        if read.reference_end + INTPRIM_SEARCH_WINDOW < genome.get_reference_length(read.reference_name):
            end_seq = genome.fetch(read.reference_name, read.reference_end, read.reference_end + INTPRIM_SEARCH_WINDOW).upper()
            right_intprim = check_intprim(end_seq)
        return left_intprim, right_intprim

    def _get_both_polyA(read):
        left_polyA, right_polyA = 0, 0
        read_seq = read.query_sequence
        if read.cigartuples[0][0] == pysam.CIGAR_OPS.CSOFT_CLIP:  # check left polyA
            end_seq = get_reverse_complement(read_seq[:read.cigartuples[0][1]])
            left_polyA = check_polyA(end_seq)
        if read.cigartuples[-1][0] == pysam.CIGAR_OPS.CSOFT_CLIP:  # check right polyA
            end_seq = read_seq[-1 * read.cigartuples[-1][1]:]
            right_polyA = check_polyA(end_seq)
        return left_polyA, right_polyA

    @classmethod
    def from_read(cls, read, junc_direction=None, *, genome=None):
        """Create a ReadRec from a pysam aligned read."""
        # FIXME switch to pycbio.hgdata.cigar
        align_start = read.reference_start
        ref_pos = align_start
        intron_blocks = []
        has_match = False
        for block in read.cigartuples:
            if block[0] == pysam.CIGAR_OPS.CREF_SKIP:  # intron
                if has_match:
                    intron_blocks.append([ref_pos, ref_pos + block[1]])
                # this fixes weird bug if there's an intron, then an insertion, then another intron???
                elif len(intron_blocks) > 0:
                    intron_blocks[-1][1] += block[1]
                has_match = False
                ref_pos += block[1]
            elif block[0] in (pysam.CIGAR_OPS.CMATCH, pysam.CIGAR_OPS.CEQUAL, pysam.CIGAR_OPS.CDIFF, pysam.CIGAR_OPS.CDEL):  # consumes reference
                ref_pos += block[1]
                if block[0] in (pysam.CIGAR_OPS.CMATCH, pysam.CIGAR_OPS.CEQUAL, pysam.CIGAR_OPS.CDIFF):
                    has_match = True
        if junc_direction not in {'+', '-'}:
            junc_direction = "-" if read.is_reverse else "+"
        juncs = tuple(Junc(blk[0], blk[1]) for blk in intron_blocks)
        left_polyA, right_polyA = cls._get_both_polyA(read)
        left_intprim, right_intprim = 0, 0
        if genome is not None:
            left_intprim, right_intprim = cls._get_both_intprim(read, genome)

        return cls(read.reference_name, junc_direction, juncs, align_start, ref_pos, read.query_name, polyA=(left_polyA, right_polyA), intprim=(left_intprim, right_intprim))

    @classmethod
    def from_junctions(cls, chrom, start, end, name, score, strand, juncs, *, polyA=None, intprim=None):
        """Create a ReadRec from junction coordinates."""
        return cls(chrom, strand, tuple(juncs), start, end, name, score=score, polyA=polyA, intprim=intprim)


class IsoWithReads:
    # keep on one copy of junction chain
    _juncs_cache = {}

    @classmethod
    def _intern_juncs(cls, juncs):
        return cls._juncs_cache.setdefault(juncs, juncs)

    def __init__(self, chrom, strand, juncs, start=None, end=None, reads=None, gene_id=None, transcript_id=None):
        self.chrom = chrom
        self.strand = strand
        self.juncs = self._intern_juncs(juncs)
        self.start = start
        self.end = end
        self._name = None
        self._score = None
        self.reads = reads if reads is not None else []
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.end5confidence = None
        self.end3confidence = None

    @property
    def name(self):
        if self._name is None:
            if self.start is None:
                return 'FLISO' + str(abs(hash(tuple(self.juncs))))
            else:
                return 'FLISO' + str(abs(hash(tuple(self.exons))))
        else:
            return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def exons(self):
        if self.start is None:
            self.start = int(median(self.starts))
            self.end = int(median(self.ends))
        return get_exons(self)

    @property
    def starts(self):
        return [x.start for x in self.reads]

    @property
    def ends(self):
        return [x.end for x in self.reads]

    @property
    def num_reads(self):
        return len(self.reads)

    @property
    def score(self):
        if self._score is None:
            return len(self.reads)
        else:
            return self._score

    @score.setter
    def score(self, new_score):
        self._score = new_score

    @property
    def genomic_length(self):
        if self.start is None or self.end is None:
            return None
        return self.end - self.start

    def reset_from_exons(self, exons):
        """Update ReadRec from a list of Exon objects."""
        self.start = exons[0].start
        self.end = exons[-1].end
        self.juncs = tuple(exons_to_juncs(sorted(exons)))

    def get_sequence(self, genome):
        if self.start is None or self.end is None:
            return None
        return get_sequence_for_exons(genome, self.chrom, self.strand, self.exons)

    @classmethod
    def from_readrec(cls, readrec):
        return cls(readrec.chrom, readrec.strand, readrec.juncs)

    @classmethod
    def from_other(cls, iso, newstart=None, newend=None, newreads=[], newstrand=None):
        if newstrand is None:
            newstrand = iso.strand
        return cls(iso.chrom, newstrand, iso.juncs, newstart, newend, newreads, iso.gene_id, iso.transcript_id)
