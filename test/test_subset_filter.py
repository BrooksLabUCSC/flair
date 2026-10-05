from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import filter_spliced_iso
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
                                 junc_to_names, isos, {}, '+', end_window=100)
    return kept


def test_subset_of_one_isoform_containing_both_terminal_exons():
    # an extra exon on each side, the candidate's terminal exons inside internal exons
    both = _iso('both', (Junc(900, 950),) + JUNCS + (Junc(1750, 1800),), 800, 1900, 50)
    assert not _is_kept([both])


def test_not_a_subset_of_two_isoforms_each_containing_one_side():
    # one shares the first exon's splice site and extends past the last exon, the other
    # the reverse; neither contains the candidate
    left = _iso('left', (Junc(900, 950),) + JUNCS + (Junc(1600, 1650),), 800, 1800, 50)
    right = _iso('right', (Junc(1100, 1150),) + JUNCS + (Junc(1750, 1800),), 1050, 1900, 50)
    assert _is_kept([left, right])


def _annots(*transcripts):
    "annotation with transcripts as (start, end, tags) on +"
    annots = AnnotData()
    for start, end, tags in transcripts:
        _save_transcript_ends(annots, '+', start, end, tags)
    for side in (annots.transcript_ends, annots.confirmed_starts, annots.confirmed_ends):
        side['+'].sort()
    annots.junc_to_gene = {}
    return annots


def test_truncated_on_one_side_needs_an_annotated_end_only_there():
    # shares the superset's first exon, so truncated only at its last exon (1700)
    superset = _iso('sup', JUNCS + (Junc(1750, 1800),), 1000, 1900, 50)
    assert not _is_kept([superset], _annots())
    assert _is_kept([superset], _annots((500, 1690, [])))  # an annotated end at 1690
    assert not _is_kept([superset], _annots((990, 3000, [])))  # only an annotated start matches
    # an end that wasn't found doesn't count
    assert not _is_kept([superset], _annots((500, 1690, ['mRNA_end_NF'])))


def test_truncated_on_both_sides_needs_one_transcripts_two_ends():
    both = _iso('both', (Junc(900, 950),) + JUNCS + (Junc(1750, 1800),), 800, 1900, 50)
    assert not _is_kept([both], _annots((990, 3000, []), (500, 1690, [])))
    assert _is_kept([both], _annots((990, 1690, [])))
