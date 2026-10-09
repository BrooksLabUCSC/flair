from types import SimpleNamespace
from flair.isoform_data import Exon, Isoform, Junc
from flair.flair_transcriptome import filter_spliced_iso, promote_backup_subsets
from flair.annotation_data import AnnotData, _save_transcript_ends

# the candidate: exons 1000-1200, 1300-1400 and 1500-1700
JUNCS = (Junc(1200, 1300), Junc(1400, 1500))
ANNOTS = SimpleNamespace(junc_to_gene={}, has_transcript_ends=lambda *args: False,
                         has_transcript_end=lambda *args: False)


def _iso(name, juncs, start, end, reads):
    iso = Isoform('chr1', '+', juncs, start, end)
    iso.name = name
    iso.reads.extend([None] * reads)
    return iso


def _is_kept(others, annots=ANNOTS):
    cand = _iso('cand', JUNCS, 1000, 1700, 10)
    isos = {iso.name: iso for iso in [cand] + others}
    junc_to_names = {}
    for iso in isos.values():
        for j in iso.juncs:
            junc_to_names.setdefault(j, set()).add(iso.name)
    kept, _ = filter_spliced_iso('nosubset', 1, cand.juncs, cand.exons, 'cand', cand.num_reads, annots,
                                 junc_to_names, isos, {}, '+', end_window=100, ss_window=10)
    return kept


def test_subset_of_one_isoform_containing_both_terminal_exons():
    # an extra exon on each side, the candidate's terminal exons inside internal exons
    both = _iso('both', (Junc(900, 950),) + JUNCS + (Junc(1750, 1800),), 800, 1900, 50)
    assert not _is_kept([both])


def test_the_supersets_making_a_subset_are_recorded():
    both = _iso('both', (Junc(900, 950),) + JUNCS + (Junc(1750, 1800),), 800, 1900, 50)
    other = _iso('other', (Junc(100, 200),), 50, 300, 50)
    isos = {iso.name: iso for iso in [_iso('cand', JUNCS, 1000, 1700, 10), both, other]}
    junc_to_names = {}
    for iso in isos.values():
        for j in iso.juncs:
            junc_to_names.setdefault(j, set()).add(iso.name)
    supersets = []
    kept, _ = filter_spliced_iso('nosubset', 1, JUNCS, isos['cand'].exons, 'cand', 10, ANNOTS, junc_to_names, isos, {}, '+',
                                 end_window=100, ss_window=10, supersets=supersets)
    assert not kept and supersets == ['both']


def test_not_a_subset_of_two_isoforms_each_containing_one_side():
    # one shares the first exon's splice site and extends past the last exon, the other
    # the reverse; neither contains the candidate
    left = _iso('left', (Junc(900, 950),) + JUNCS + (Junc(1600, 1650),), 800, 1800, 50)
    right = _iso('right', (Junc(1100, 1150),) + JUNCS + (Junc(1750, 1800),), 1050, 1900, 50)
    assert _is_kept([left, right])


def _annots(*transcripts):
    "annotation with + strand transcripts given as (exons, tags)"
    annots = AnnotData()
    for exons, tags in transcripts:
        exons = [Exon(s, e) for s, e in exons]
        _save_transcript_ends(annots, '+', exons[0].start, exons[-1].end, tags, exons)
    annots.transcript_ends['+'].sort()
    for ends in annots.confirmed_terminal_ends.values():
        ends.sort()
    annots.junc_to_gene = {}
    return annots


def test_truncated_on_one_side_needs_an_annotated_end_there_with_the_same_exon():
    # shares the superset's first exon, so truncated only at its last exon (1500-1700)
    superset = _iso('sup', JUNCS + (Junc(1750, 1800),), 1000, 1900, 50)
    assert not _is_kept([superset], _annots())
    # an annotated last exon 1500-1690
    assert _is_kept([superset], _annots(([(500, 1200), (1500, 1690)], ['basic'])))
    # a non-basic transcript's end, as a retained_intron fragment's, isn't enough
    assert not _is_kept([superset], _annots(([(500, 1200), (1500, 1690)], [])))
    # but one non-basic transcript matching both ends still is
    assert _is_kept([superset], _annots(([(990, 1200), (1500, 1690)], [])))
    # a last exon starting within ss_window of the subset's, 1506, and not beyond it
    assert _is_kept([superset], _annots(([(500, 1200), (1506, 1690)], ['basic'])))
    assert not _is_kept([superset], _annots(([(500, 1200), (1520, 1690)], ['basic'])))
    # an annotated end at 1690, but of a last exon starting elsewhere, as a fragment's
    assert not _is_kept([superset], _annots(([(500, 1450), (1600, 1690)], ['basic'])))
    # only an annotated start matches
    assert not _is_kept([superset], _annots(([(990, 1200), (1500, 3000)], [])))
    # an end that wasn't found doesn't count
    assert not _is_kept([superset], _annots(([(500, 1200), (1500, 1690)], ['basic', 'mRNA_end_NF'])))


def test_truncated_on_both_sides_needs_one_transcripts_two_ends():
    both = _iso('both', (Junc(900, 950),) + JUNCS + (Junc(1750, 1800),), 800, 1900, 50)
    assert not _is_kept([both], _annots(([(990, 1200), (1500, 3000)], []), ([(500, 1200), (1500, 1690)], [])))
    assert _is_kept([both], _annots(([(990, 1200), (1500, 1690)], [])))


###
# backup subsets: a subset reported after all when every isoform it is a subset of fails support
###
BACKUP_ARGS = SimpleNamespace(subset_backup_support=10, sjc_support=1, single_exon_support=3, frac_support=0.05)


def _read(name, start, end):
    return SimpleNamespace(name=name, start=start, end=end)


def _backup_case(superset_reads, subset_reads=12, assigned=()):
    """a kept isoform 'big' (100 reads) of gene g, the superset 'sup' with
    superset_reads reads, and the removed subset with subset_reads reads; the
    ends file lists the assigned reads"""
    big = _iso('big', (Junc(5000, 5100),), 4900, 5300, 0)
    sup = _iso('sup', (Junc(900, 950),) + JUNCS, 800, 1700, 0)
    for iso in (big, sup):
        iso.gene_id, iso.unpadded_ends = 'g', (iso.start, iso.end)
    subset = _iso('cand', JUNCS, 1000, 1700, 0)
    subset.reads[:] = [_read(f'r{i}', 1000 + i, 1700 - i) for i in range(subset_reads)]
    final = {'big': big, 'sup': sup}
    counts = {'big': [100, 100], 'sup': [superset_reads, superset_reads]}
    gene_to_tot = {'g': [100 + superset_reads] * 3}
    return final, counts, gene_to_tot, {'cand': (subset, [sup])}


def _promote(tmp_path, final, counts, gene_to_tot, removed, assigned=()):
    ends = tmp_path / 'ends.tsv'
    ends.write_text(''.join(f'{name}\tsup\n' for name in assigned))
    promote_backup_subsets(BACKUP_ARGS, removed, final, counts, gene_to_tot, str(ends), None)
    return [iso for iso in final.values() if iso.backup_subset]


def test_subset_promoted_when_its_superset_fails_support(tmp_path):
    # the superset's 2 of 102 reads fail frac_support 0.05
    final, counts, gene_to_tot, removed = _backup_case(2)
    (promoted,) = _promote(tmp_path, final, counts, gene_to_tot, removed)
    assert promoted.juncs == JUNCS and promoted.gene_id == 'g' and counts[promoted.name] == [12, 12]
    assert gene_to_tot['g'] == [114] * 3 and promoted.unpadded_ends == (promoted.start, promoted.end)


def test_subset_not_promoted_when_its_superset_passes(tmp_path):
    assert _promote(tmp_path, *_backup_case(30)) == []


def test_subset_promoted_only_with_subset_backup_support_unassigned_reads(tmp_path):
    # 12 reads, 3 of them assigned to another isoform: 9 left, under 10
    final, counts, gene_to_tot, removed = _backup_case(2)
    assert _promote(tmp_path, final, counts, gene_to_tot, removed, assigned=('r0', 'r1', 'r2')) == []


def test_a_fragment_of_a_promoted_subset_is_not_promoted(tmp_path):
    # 'frag' is a subset of 'cand', itself promoted; listed first, it is still decided after
    final, counts, gene_to_tot, removed = _backup_case(2)
    cand, sup = removed['cand'][0], removed['cand'][1][0]
    frag = _iso('frag', JUNCS[1:], 1350, 1700, 0)
    frag.reads[:] = [_read(f'f{i}', 1350, 1700) for i in range(12)]
    removed = {'frag': (frag, [cand, sup]), 'cand': removed['cand']}
    (promoted,) = _promote(tmp_path, final, counts, gene_to_tot, removed)
    assert promoted.juncs == JUNCS
