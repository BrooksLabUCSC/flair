import random
from types import SimpleNamespace
import pysam
from flair.isoform_data import Exon
from flair.read_correction import _correct_and_group_read

# _correct_and_group_read, for a read whose annotation match lacks one of its introns

HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 10000}]})
# an annotated transcript with exons 1000-1100 and 1300-1700: one junction, 1100-1300
ANNOTS = SimpleNamespace(transcript_to_exons={('TX', 'GENE'): (Exon(1000, 1100), Exon(1300, 1700))},
                         gene_to_strand={'GENE': '+'})


def _read(cigar):
    read = pysam.AlignedSegment(HEADER)
    read.query_name, read.reference_id, read.reference_start = 'r', 0, 1000
    read.cigarstring = cigar
    read.query_sequence = 'ACGT' * (read.query_alignment_length // 4) + 'ACGT'[:read.query_alignment_length % 4]
    return read


class _Corrector:
    "junction correction that finds the read's own junctions supported, or not"
    def __init__(self, supported):
        self.supported, self.called = supported, False

    def correct_readrec(self, readrec, trust_strand, genome):
        self.called = True
        return self.supported


def _correct(read, corrector):
    sj_to_ends = {}
    # matched to TX from its first junction to its first, 100 bases before it and 400 past
    _correct_and_group_read(read, read_to_annot_transcript={'r': ('TX_GENE', 0, 100, 0, 400)}, annots=ANNOTS,
                            junction_corrector=corrector, sj_to_ends=sj_to_ends, genome=None,
                            keep_single_exon=True, trust_strand=False, check_motifs=False)
    ((chrom, juncs), isoform), = sj_to_ends.items()
    read_rec = isoform.reads[0]
    return tuple((j.start, j.end) for j in juncs), (read_rec.start, read_rec.end)


def test_a_supported_intron_inside_an_annotated_exon_is_kept():
    # spliced at 1100-1300 and, inside the annotated last exon, at 1400-1600
    corrector = _Corrector(supported=True)
    assert _correct(_read('100M200N100M200N100M'), corrector)[0] == ((1100, 1300), (1400, 1600))
    assert corrector.called


def test_an_unsupported_extra_junction_takes_the_annotation():
    juncs, ends = _correct(_read('100M200N100M200N100M'), _Corrector(supported=False))
    assert juncs == ((1100, 1300),) and ends == (1000, 1700)


def test_a_match_losing_no_junction_is_used_without_correction():
    corrector = _Corrector(supported=True)
    assert _correct(_read('100M200N400M'), corrector)[0] == ((1100, 1300),)
    assert not corrector.called


# the junctions_moved and clean_splice_sites flags correction sets on a read
random.seed(1)
SEQ = ''.join(random.choice('ACGT') for _ in range(3000))
GENOME = SimpleNamespace(fetch=lambda chrom, start, end: SEQ[start:end], get_reference_length=lambda chrom: len(SEQ))


def _genome_read(cigar, mismatches=()):
    "a read at 1000 aligned with the genome's sequence, but for mismatches at the given positions"
    read = _read(cigar)
    pos, seq = 1000, ''
    for op, length in read.cigartuples:
        if op == pysam.CMATCH:
            seq += SEQ[pos:pos + length]
        pos += length
    for m in mismatches:
        i = m - 1000 - (200 if m > 1300 else 0)
        seq = seq[:i] + ('A' if seq[i] != 'A' else 'C') + seq[i + 1:]
    read.query_sequence = seq
    return read


def _flags(read, annotated):
    sj_to_ends = {}
    _correct_and_group_read(read, read_to_annot_transcript={'r': ('TX_GENE', 0, 100, 0, 400)} if annotated else {},
                            annots=ANNOTS, junction_corrector=_Corrector(supported=True), sj_to_ends=sj_to_ends,
                            genome=GENOME, keep_single_exon=True, trust_strand=False, check_motifs=False)
    read_rec = next(iter(sj_to_ends.values())).reads[0]
    return read_rec.junctions_moved, read_rec.clean_splice_sites


def test_junctions_kept_from_the_alignment_are_checked_for_clean_splice_sites():
    # corrected from intron support without changing the junction
    assert _flags(_genome_read('100M200N400M'), annotated=False) == (False, True)
    # two mismatches next to the splice site
    assert _flags(_genome_read('100M200N400M', mismatches=(1098, 1099)), annotated=False) == (False, False)


def test_junctions_from_an_annotation_match_are_moved_when_they_differ():
    # the annotation's junction is the read's own: not moved
    assert _flags(_genome_read('100M200N400M'), annotated=True) == (False, True)
    # a read missing the annotated junction is given it: moved, and not checked
    assert _flags(_genome_read('700M'), annotated=True) == (True, None)
