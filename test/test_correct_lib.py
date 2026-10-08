import pytest
from flair.pycbio.hgdata.bed import Bed
from flair.isoform_data import ReadRec, Junc
from flair.intron_support import IntronSupport
from flair.junction_correct import JunctionCorrector, UNKNOWN_STRAND

# lots of long lines for BEDs
# flake8: noqa: E501

###
# tests for loading intron support data
###
def _assert_introns(over, expect):
    # order-insensitive: overlap() has no defined order
    assert sorted(str(o) for o in over) == sorted(expect)

def _basic_load_reads_annot_support_test(intron_support):
    "same data, loaded from BED or STAR SJ"
    # chr12:25209911-25215436(-)
    start = 25209911
    end = 25215436
    over = intron_support.overlap("chr12", start, end, 5)
    _assert_introns(over, [
        'SupportIntron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
    ])

def test_load_bed_annot_support(basic_gtf_data):
    intron_support = IntronSupport()
    intron_support.load_introns_bed("input/basic.shortread_junctions.bed")
    intron_support.load_gtf(basic_gtf_data)
    _basic_load_reads_annot_support_test(intron_support)

def test_load_star_annot_support(basic_gtf_data):
    intron_support = IntronSupport()
    intron_support.load_star("input/basic.shortread_junctions.tab")
    intron_support.load_gtf(basic_gtf_data)
    _basic_load_reads_annot_support_test(intron_support)


###
# tests for correcting introns
###
@pytest.fixture(scope="session")
def basic_support(basic_gtf_data):
    intron_support = IntronSupport()
    intron_support.load_introns_bed("input/basic.shortread_junctions.bed")
    intron_support.load_gtf(basic_gtf_data)
    return intron_support

@pytest.fixture(scope="session")
def basic_corrector(basic_support):
    # standard parameters for flair_correct
    return JunctionCorrector(basic_support, 15, 1)

def _readrec_from_bed_str(bed_str):
    """Parse a BED string into a ReadRec for testing."""
    bed = Bed.parse(bed_str.split('\t'))
    juncs = tuple(Junc(bed.blocks[i].end, bed.blocks[i + 1].start)
                  for i in range(len(bed.blocks) - 1))
    return ReadRec(bed.chrom, bed.strand, juncs, bed.chromStart, bed.chromEnd, bed.name)

def test_no_support(basic_corrector):
    readrec = _readrec_from_bed_str("chr20	35542150	35557051	HISEQ:1287:HKCG7BCX3:1:1106:16788:27192	60	+	35542150	35557051	27,158,119	10	35,71,88,120,94,166,62,132,65,79,	0,172,362,671,5261,6358,6657,13847,14056,14822,")
    assert basic_corrector.correct_readrec(readrec) is False

def test_adjust_start(basic_corrector):
    # start of intron adjusted
    readrec = _readrec_from_bed_str("chr20	35548585	35557571	HISEQ:1287:HKCG7BCX3:1:1106:18816:99230	60	+	35548585	35557571	27,158,119	7	89,58,32,97,65,81,132,	0,222,6458,7447,7621,8387,8854,")
    expt = _readrec_from_bed_str("chr20	35548585	35557571	HISEQ:1287:HKCG7BCX3:1:1106:18816:99230	60	+	35548585	35557571	27,158,119	7	89,58,32,97,65,83,132,	0,222,6458,7447,7621,8387,8854,")
    assert basic_corrector.correct_readrec(readrec) is True
    assert readrec.juncs == expt.juncs
    assert readrec.strand == expt.strand

def test_adjust_end(basic_corrector):
    # end of intron adjusted
    readrec = _readrec_from_bed_str("chr20	35548508	35557489	HISEQ:1287:HKCG7BCX3:1:1107:10927:56503	60	+	35548508	35557489	27,158,119	7	166,58,32,97,65,87,46,	0,299,6535,7524,7698,8464,8935,")
    expt = _readrec_from_bed_str("chr20	35548508	35557489	HISEQ:1287:HKCG7BCX3:1:1107:10927:56503	60	+	35548508	35557489	27,158,119	7	166,58,32,97,65,83,50,	0,299,6535,7524,7698,8464,8931,")
    assert basic_corrector.correct_readrec(readrec) is True
    assert readrec.juncs == expt.juncs
    assert readrec.strand == expt.strand

def test_adjust_both(basic_corrector):
    # both of intron adjusted
    readrec = _readrec_from_bed_str("chr17	64499205	64504277	HISEQ:1287:HKCG7BCX3:1:1101:17977:51378	60	-	64499205	64504277	217,95,2	11	1121,225,64,62,111,173,161,142,66,134,56,	0,1343,2800,2956,3233,3720,3982,4224,4597,4777,5016,")
    expt = _readrec_from_bed_str("chr17	64499205	64504277	HISEQ:1287:HKCG7BCX3:1:1101:17977:51378	60	-	64499205	64504277	217,95,2	11	1121,225,60,62,111,173,161,142,66,134,56,	0,1343,2804,2956,3233,3720,3982,4224,4597,4777,5016,")
    assert basic_corrector.correct_readrec(readrec) is True
    assert readrec.juncs == expt.juncs
    assert readrec.strand == expt.strand

###
# trust_strand: the read's strand is kept and only introns on it are used
###
def _stranded_corrector(*supports):
    intron_support = IntronSupport()
    for start, end, strand in supports:
        intron_support.add_support("chr1", start, end, strand, 10)
    return JunctionCorrector(intron_support, 15, 1)

def _plus_read():
    return ReadRec("chr1", '+', (Junc(1000, 2000),), 500, 2500, "read1")

def test_opposite_strand_support():
    # like an intron-prospector junction with a semi-canonical motif on the other strand
    corrector = _stranded_corrector((1000, 2000, '-'))
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec) is True
    assert readrec.strand == '-'
    assert corrector.correct_readrec(_plus_read(), trust_strand=True) is False

def test_trust_strand_uses_same_strand_support():
    # the closest support is on the other strand, so it is only used without trust_strand
    corrector = _stranded_corrector((1000, 2000, '-'), (1004, 2004, '+'))
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec) is True
    assert (readrec.juncs, readrec.strand) == ((Junc(1000, 2000),), '-')
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec, trust_strand=True) is True
    assert (readrec.juncs, readrec.strand) == ((Junc(1004, 2004),), '+')

def test_trust_strand_unknown_strand_support():
    # support of unknown strand can't give a strand, but the read's can
    corrector = _stranded_corrector((1000, 2000, '.'))
    assert corrector.correct_readrec(_plus_read()) is False
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec, trust_strand=True) is True
    assert readrec.strand == '+'

###
# junction motifs as strand evidence: without an annotation or a canonical GT-AG
# motif, an intron's strand is weak evidence and the read's strand is left unknown
###
class _FakeGenome:
    "a chr1 of N, with the given sequence placed at positions"
    def __init__(self, placed):
        seq = ['N'] * 3000
        for pos, bases in placed.items():
            seq[pos:pos + len(bases)] = bases
        self.seq = ''.join(seq)

    def fetch(self, chrom, start, end):
        return self.seq[start:end]

def _motif_genome(donor, acceptor):
    # the intron of _plus_read is 1000-2000
    return _FakeGenome({1000: donor, 1998: acceptor})

def test_canonical_motif_gives_strand():
    corrector = _stranded_corrector((1000, 2000, '-'))
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec, genome=_motif_genome('CT', 'AC')) is True
    assert readrec.strand == '-'

def test_semicanonical_motif_gives_unknown_strand():
    # GC-AG on the - strand, as intron-prospector reports it
    corrector = _stranded_corrector((1000, 2000, '-'))
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec, genome=_motif_genome('CT', 'GC')) is True
    assert readrec.strand == UNKNOWN_STRAND
    # without the genome the motif isn't checked
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec) is True
    assert readrec.strand == '-'

def test_annotated_semicanonical_intron_gives_strand():
    intron_support = IntronSupport()
    intron_support.add_support("chr1", 1000, 2000, '-')
    corrector = JunctionCorrector(intron_support, 15, 1)
    readrec = _plus_read()
    assert corrector.correct_readrec(readrec, genome=_motif_genome('CT', 'GC')) is True
    assert readrec.strand == '-'

def test_one_canonical_intron_strands_the_read():
    corrector = _stranded_corrector((1000, 2000, '-'), (2100, 2400, '-'))
    readrec = ReadRec("chr1", '+', (Junc(1000, 2000), Junc(2100, 2400)), 500, 2900, "read2")
    genome = _FakeGenome({1000: 'CT', 1998: 'GC', 2100: 'CT', 2398: 'AC'})
    assert corrector.correct_readrec(readrec, genome=genome) is True
    assert readrec.strand == '-'

def _strong_plus_weak_minus_read():
    return ReadRec("chr1", '+', (Junc(1000, 2000), Junc(2100, 2400)), 500, 2900, "read3")

# GT-AG on + at 1000-2000; GC-AG on - at 2100-2400; GT-AG on + at 2104-2404
_MIXED_GENOME = _FakeGenome({1000: 'GT', 1998: 'AG', 2100: 'CT', 2398: 'GC', 2104: 'GT', 2402: 'AG'})

def test_weak_opposite_strand_recorrected_on_strong_strand():
    corrector = _stranded_corrector((1000, 2000, '+'), (2100, 2400, '-'), (2104, 2404, '+'))
    readrec = _strong_plus_weak_minus_read()
    assert corrector.correct_readrec(readrec, genome=_MIXED_GENOME) is True
    assert (readrec.juncs, readrec.strand) == ((Junc(1000, 2000), Junc(2104, 2404)), '+')
    # without the motif check, the strands conflict
    assert corrector.correct_readrec(_strong_plus_weak_minus_read()) is False

def test_weak_opposite_strand_without_strong_strand_support():
    corrector = _stranded_corrector((1000, 2000, '+'), (2100, 2400, '-'))
    assert corrector.correct_readrec(_strong_plus_weak_minus_read(), genome=_MIXED_GENOME) is False
