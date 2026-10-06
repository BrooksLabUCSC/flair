import io
from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.direct_assignment import covers_unique_bounds, direct_assignments, write_direct_assignments

# exons 1000-1100, 1300-1400 and 1600-1700
JUNCS = (Junc(1100, 1300), Junc(1400, 1600))


def _isoform(name, juncs=JUNCS, strand='+'):
    iso = Isoform('chr1', strand, juncs, 900, 1800)
    iso.name = name
    return iso


def _read(name, juncs=JUNCS, start=1000, end=1700, moved=False, clean=True, annotation=False, terminal_moved=False):
    return SimpleNamespace(name=name, juncs=juncs, start=start, end=end, junctions_moved=moved, clean_splice_sites=clean,
                           junctions_from_annotation=annotation, terminal_junctions_moved=terminal_moved)


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
