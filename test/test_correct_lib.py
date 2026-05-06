import pytest
from flair.pycbio.hgdata.bed import Bed
from flair.isoform_data import ReadRec, Junc
from flair.intron_support import IntronSupport
from flair.junction_correct import JunctionCorrector

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
    over = intron_support.overlap("chr12", start, start + 1)
    _assert_introns(over, [
        'SupportIntron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
        'SupportIntron(chr12:25209911-25225613(-) annot=True read=False read_cnt=0',
        'SupportIntron(chr12:25209911-25245273(-) annot=True read=False read_cnt=0',
        'SupportIntron(chr12:25209911-25225613(-) annot=False read=True read_cnt=96',
    ])

    over = intron_support.overlap("chr12", end - 1, end)
    _assert_introns(over, [
        'SupportIntron(chr12:25209911-25215436(-) annot=True read=True read_cnt=12',
    ])

    over = intron_support.overlap_introns("chr12", start, end, 5)
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
