from types import SimpleNamespace
from flair.isoform_data import Junc, Exon
from flair.flair_transcriptome import _group_novel_by_splicing, _group_by_read_span


def _iso(start, end, juncs=()):
    juncs = tuple(Junc(s, e) for s, e in juncs)
    bounds = [start] + [x for j in juncs for x in (j.start, j.end)] + [end]
    exons = [Exon(bounds[i], bounds[i + 1]) for i in range(0, len(bounds), 2)]
    return SimpleNamespace(start=start, end=end, juncs=juncs, exons=exons)


def _groups(group_func, isos):
    return sorted(sorted(g) for g in group_func(list(isos), isos))


# two spliced loci, the first with a long read tail reaching into the second
LOCI = {'a': _iso(1000, 9000, [(1200, 1800)]),      # tail from 1900 to 9000
        'b': _iso(8000, 9500, [(8200, 8800)])}


def test_read_tail_bridges_loci_only_with_read_spans():
    assert _groups(_group_by_read_span, LOCI) == [['a', 'b']]
    assert _groups(_group_novel_by_splicing, LOCI) == [['a'], ['b']]


def test_shared_splice_site_groups():
    # c and d share only the 3' splice site at 1800, and their junction spans overlap anyway;
    # e shares nothing and its junction span is elsewhere
    isos = {'c': _iso(1000, 2000, [(1200, 1800)]), 'd': _iso(1300, 2100, [(1500, 1800)]),
            'e': _iso(5000, 6000, [(5200, 5800)])}
    assert _groups(_group_novel_by_splicing, isos) == [['c', 'd'], ['e']]


def test_overlapping_junction_spans_merge_without_shared_sites():
    # f's intron contains g's; they share no splice site
    isos = {'f': _iso(1000, 5000, [(1100, 4900)]), 'g': _iso(2000, 3000, [(2100, 2900)])}
    assert _groups(_group_novel_by_splicing, isos) == [['f', 'g']]


def test_single_exon_attaches_by_exon_overlap_only():
    # s1 overlaps h's first exon; s2 lies in h's intron, overlapping no exon
    isos = {'h': _iso(1000, 5000, [(1200, 4800)]),
            's1': _iso(900, 1100), 's2': _iso(2000, 3000)}
    assert _groups(_group_novel_by_splicing, isos) == [['h', 's1'], ['s2']]
    # with read spans, s2 is inside h's span and joins it
    assert _groups(_group_by_read_span, isos) == [['h', 's1', 's2']]
