import io
from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
import pysam
from flair.isoform_data import ReadRec
from flair.direct_assignment import covers_unique_bounds, direct_assignments, unassignable_reads, write_direct_assignments

# exons 1000-1100, 1300-1400 and 1600-1700
JUNCS = (Junc(1100, 1300), Junc(1400, 1600))


def _isoform(name, juncs=JUNCS, strand='+'):
    iso = Isoform('chr1', strand, juncs, 900, 1800)
    iso.name = name
    return iso


def _read(name, juncs=JUNCS, start=1000, end=1700, moved=False, clean=True, annotation=False, terminal_moved=False,
          clipping=(0, 0)):
    return SimpleNamespace(name=name, juncs=juncs, start=start, end=end, junctions_moved=moved, clean_splice_sites=clean,
                           junctions_from_annotation=annotation, terminal_junctions_moved=terminal_moved, clipping=clipping)


def _assigned(reads, isoforms, bounds=None):
    return {n: i.name for n, (i, _) in direct_assignments(reads, isoforms, bounds or {}).items()}


def test_an_exact_clean_match_is_assigned_directly():
    assert _assigned([_read('r')], [_isoform('tx')]) == {'r': 'tx'}


def test_moved_or_unclean_junctions_are_realigned():
    assert _assigned([_read('moved', moved=True), _read('unclean', clean=False), _read('unchecked', clean=None)],
                     [_isoform('tx')]) == {}


def test_annotation_corrected_reads_need_their_terminal_junctions_unmoved():
    # a missed short exon added inside the read: assigned, whatever its splice sites
    internal = _read('internal', moved=True, clean=None, annotation=True)
    # a junction added or dropped at an end of the read
    terminal = _read('terminal', moved=True, clean=None, annotation=True, terminal_moved=True)
    # moved by intron support rather than an annotation match
    introns = _read('introns', moved=True, clean=None)
    assert _assigned([internal, terminal, introns], [_isoform('tx')]) == {'internal': 'tx'}


def test_a_chain_of_two_isoforms_is_realigned():
    assert _assigned([_read('r')], [_isoform('tx1'), _isoform('tx2')]) == {}


def test_unique_bounds_must_be_covered():
    plus, minus = _isoform('tx'), _isoform('tx', strand='-')
    # on +, side 0 is the genomic start: the read must start 50 + 5 bases before the first junction
    assert covers_unique_bounds(plus, {0: 50}, 1045, 1700)
    assert not covers_unique_bounds(plus, {0: 50}, 1050, 1700)
    # on -, side 0 is the genomic end, past the last junction
    assert covers_unique_bounds(minus, {0: 50}, 1000, 1655)
    assert not covers_unique_bounds(minus, {0: 50}, 1000, 1650)
    assert _assigned([_read('r', start=1060)], [plus], {'tx': {0: 50}}) == {}


def test_ends_are_written_as_count_sam_transcripts_writes_them():
    fh = io.StringIO()
    write_direct_assignments({'r': (_isoform('tx'), _read('r', start=1010, end=1690))}, fh)
    assert fh.getvalue() == 'r\ttx\t0\t90\t110\t1\t90\t110\n'


def test_reads_no_kept_isoform_can_take_are_not_realigned():
    isoforms = [_isoform('tx')]
    # its chain is the kept isoform's
    kept = _read('kept')
    # a chain no kept isoform has
    other = _read('other', juncs=(Junc(1100, 1250),), end=1400)
    # the kept isoform's first junction only, ending 2 bases short of its next
    # splice site (1400) and clipped there, as when a short last exon wasn't aligned
    clipped = _read('clipped', juncs=(Junc(1100, 1300),), end=1398, clipping=(0, 12))
    # the same, but not clipped, or ending too far from the splice site
    unclipped = _read('unclipped', juncs=(Junc(1100, 1300),), end=1398)
    too_far = _read('too_far', juncs=(Junc(1100, 1300),), end=1350, clipping=(0, 12))
    reads = [kept, other, clipped, unclipped, too_far]
    assert unassignable_reads(reads, isoforms, 10) == {'other', 'unclipped', 'too_far'}


def test_a_read_records_its_alignments_clipping():
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chr1', 'LN': 10000}]})
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id, read.reference_start = 'r', 0, 1000
    read.cigarstring, read.query_sequence = '2S100M200N100M3S', 'A' * 205
    assert ReadRec.from_read(read).clipping == (2, 3)
