"""
Annotation and orthogonal intron support.
"""
import sys
from collections import defaultdict
from flair import MIN_INTRON_SIZE, MAX_INTRON_SIZE, FlairInputDataError
from flair.gtf_io import GtfData
from flair.interval_index import IntervalIndex
from flair.pycbio.hgdata.bed import BedReader
from flair.pycbio.tsv import TsvReader

_VALID_STRANDS = {'+', '-', '.'}

class SupportIntron:
    """Coordinates and support.  Strand could None if not available, but currently
    enforce having strand."""
    __slots__ = ("chrom", "start", "end", "strand",
                 "annot_supported", "read_supported", "read_support_cnt")

    def __init__(self, chrom, start, end, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.annot_supported = self.read_supported = False
        self.read_support_cnt = 0

    def __str__(self):
        return (f"SupportIntron({self.chrom}:{self.start}-{self.end}({self.strand}) "
                f"annot={self.annot_supported} read={self.read_supported} read_cnt={self.read_support_cnt}")

class IntronSupport:
    """
    Table of intron support index by both start and end positions
    """
    def __init__(self, *, min_intron_size=MIN_INTRON_SIZE, max_intron_size=MAX_INTRON_SIZE):
        # dict index by chrom of per-chrom interval indexes, keyed on first base of donor and last base of the acceptor sites
        self.coords_maps = defaultdict(IntervalIndex)
        self.min_intron_size = min_intron_size
        self.max_intron_size = max_intron_size
        self.chroms = set()

    def _find_point_strand(self, chrom, point, strand):
        for intron in self.coords_maps[chrom].overlap(point, point + 1):
            if intron.strand == strand:
                return intron
        return None

    def _find_intron(self, chrom, start, end, strand):
        donor = self._find_point_strand(chrom, start, strand)
        if donor is not None:
            accept = self._find_point_strand(chrom, end - 1, strand)
            if accept is donor:
                return donor  # same intron
        return None

    def _add_intron(self, chrom, start, end, strand):
        assert start < end
        intron = SupportIntron(chrom, start, end, strand)
        self.coords_maps[chrom].add(start, start + 1, intron)
        self.coords_maps[chrom].add(end - 1, end, intron)
        self.chroms.add(chrom)
        return intron

    def _add_support(self, chrom, start, end, strand, read_count):
        intron = self._find_intron(chrom, start, end, strand)
        if intron is None:
            intron = self._add_intron(chrom, start, end, strand)
        if read_count is None:
            intron.annot_supported = True
        else:
            intron.read_supported = True
            intron.read_support_cnt += max(read_count, 1)

    def add_support(self, chrom, start, end, strand, read_count=None):
        """A read_count None indicates annot_support.  Drop intron outside of
        configured size range"""
        if self.min_intron_size <= (end - start) <= self.max_intron_size:
            self._add_support(chrom, start, end, strand, read_count)
            return True
        return False

    def overlap(self, chrom, start, end, flank_window=0):
        """Get list of overlapping introns where either ends overlaps this range with
        a +/-bp window"""
        return self.coords_maps[chrom].overlap(start, end, slack=flank_window)

    def overlap_introns(self, chrom, start, end, flank_window=0):
        """get introns were splice junctions overlap each end of this range,
        with a +/-bp window on either of ends of the range. Empty list if no hits"""
        # only return introns that hit both ends

        starts = self.overlap(chrom, start, start + 1, flank_window)
        ends = self.overlap(chrom, end - 1, end, flank_window)

        ends_ids = set([id(e) for e in ends])
        return [intron for intron in starts
                if id(intron) in ends_ids]

    def chroms(self):
        "generator for chroms"
        return self.coords_maps.keys()

    def entries(self, chrom=None):
        "generator for (chrom, start, end, intron), optionally on a chrom (introns have two entries)"
        chroms = [chrom] if chrom is not None else self.chroms()
        for chrom in chroms:
            for start, end, intron in self.coords_maps[chrom].items():
                yield chrom, start, end, intron

    def introns(self, chrom=None):
        "generator for introns, optionally on a chrom"
        seen = set()  # introns are in twice
        for _, _, _, intron in self.entries(chrom):
            if id(intron) not in seen:
                yield intron
                seen.add(id(intron))

    def subset_for_region(self, chrom, start, end):
        """Return a new IntronSupport with introns overlapping [start, end) on chrom.
        Intron support counts and flags are copied."""
        sub = IntronSupport(min_intron_size=self.min_intron_size, max_intron_size=self.max_intron_size)
        if chrom in self.chroms:
            for intron in self.introns(chrom):
                if intron.end > start and intron.start < end:
                    new_intron = sub._add_intron(intron.chrom, intron.start, intron.end, intron.strand)
                    new_intron.annot_supported = intron.annot_supported
                    new_intron.read_supported = intron.read_supported
                    new_intron.read_support_cnt = intron.read_support_cnt
        return sub

    def dump(self, fh=sys.stderr):
        print("IntronSupport:", file=fh)
        for chrom, start, end, intron in self.entries():
            print(f"{chrom}:{start}-{end}: {intron}", file=fh)

    @staticmethod
    def _no_introns_loaded_error(file_name, file_desc, chrom_filter=None):
        msg = f"No introns loaded from {file_desc} file"
        if chrom_filter is not None:
            msg += f" for chromosome `{chrom_filter}'"
        msg += f": {file_name}"
        raise FlairInputDataError(msg)

    @staticmethod
    def _bed_strand_error(bed):
        raise FlairInputDataError(f"Invalid strand `{bed.strand}' in BED must be one of: " +
                                  ", ".join([f"'{s}'" for s in _VALID_STRANDS]))

    def _load_intron_bed(self, bed, chrom_filter):
        if (chrom_filter is not None) and (bed.chrom != chrom_filter):
            return False
        if not (6 <= bed.numStdCols <= 9):
            raise FlairInputDataError(f"intron BED must have 6 to 9 columns, found {bed.numStdCols}")
        if bed.strand not in _VALID_STRANDS:
            self._bed_strand_error(bed)
        return self.add_support(bed.chrom, bed.chromStart, bed.chromEnd, bed.strand, bed.score)

    def load_introns_bed(self, bed_file, *, chrom_filter=None):
        """load introns from a BED6, with score being the counts of reads.
        Return number of introns loaded."""
        cnt = 0
        line_num = 0
        try:
            for bed in BedReader(bed_file):
                line_num += 1
                if self._load_intron_bed(bed, chrom_filter):
                    cnt += 1
        except Exception as exc:
            raise FlairInputDataError(f"parsing intron BED failed: {bed_file} line {line_num}") from exc
        if cnt == 0:
            self._no_introns_loaded_error(bed_file, "BED", chrom_filter)
        return cnt

    def _load_annot_bed(self, bed, chrom_filter):
        if (chrom_filter is not None) and (bed.chrom != chrom_filter):
            return False
        if len(bed.blocks) == 1:
            return False
        if bed.numStdCols != 12:
            raise FlairInputDataError(f"annot BED must have 12 columns, found {bed.numStdCols}")
        if bed.strand not in _VALID_STRANDS:
            self._bed_strand_error(bed)
        for i in range(len(bed.blocks) - 1):
            self.add_support(bed.chrom, bed.blocks[i].end, bed.blocks[i + 1].start, bed.strand, None)
        return True

    def load_annot_bed(self, bed_file, *, chrom_filter=None):
        """load introns from a BED12, with score staying none to register as annot
        Return number of introns loaded."""
        cnt = 0
        line_num = 0
        try:
            for bed in BedReader(bed_file):
                line_num += 1
                if self._load_annot_bed(bed, chrom_filter):
                    cnt += 1
        except Exception as exc:
            raise FlairInputDataError(f"parsing intron BED failed: {bed_file} line {line_num}") from exc

        if cnt == 0:
            self._no_introns_loaded_error(bed_file, "BED", chrom_filter)
        return cnt

    def _load_star(self, rec, chrom_filter):
        if (chrom_filter is not None) and (rec.chrom != chrom_filter):
            return False
        strand = (None, '+', '-')[rec.strand]
        if strand is None:
            return False
        return self.add_support(rec.chrom, rec.start - 1, rec.end, strand, rec.uniq_map_cnt)

    def load_star(self, sj_file, *, chrom_filter=None):
        """load introns from STAR junction file"""
        columns = ("chrom", "start", "end", "strand", "motif", "annot", "uniq_map_cnt",
                   "multi_map_cnt", "max_overhang")
        cnt = 0
        line_num = 1  # include header
        try:
            for rec in TsvReader(sj_file, columns=columns, typeMap={"chrom": str}, defaultColType=int):
                line_num += 1
                if self._load_star(rec, chrom_filter):
                    cnt += 1
        except Exception as exc:
            raise FlairInputDataError(f"parsing STAR SJ failed: {sj_file} line {line_num}") from exc
        if cnt == 0:
            self._no_introns_loaded_error(sj_file, "STAR SJ", chrom_filter)
        return cnt

    def _load_gtf_transcript(self, transcript):
        if len(transcript.exons) == 0:
            raise FlairInputDataError(f"transcript '{transcript.transcript_id}' has no exons")
        cnt = 0
        prev_end = transcript.exons[0].end
        for exon in transcript.exons[1:]:
            if self.add_support(exon.chrom, prev_end, exon.start, exon.strand):
                cnt += 1
            prev_end = exon.end
        return cnt

    def load_gtf(self, gtf_data: GtfData):
        """Load introns from an annotation GTF. Return number of introns loaded."""
        cnt = 0
        try:
            for transcript in gtf_data.transcripts:
                cnt += self._load_gtf_transcript(transcript)
        except Exception as exc:
            raise FlairInputDataError(f"parsing annotation GTF failed: {gtf_data.gtf_file}") from exc
        if cnt == 0:
            self._no_introns_loaded_error(gtf_data.gtf_file, "GTF")
        return cnt
