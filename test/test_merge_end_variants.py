from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import filter_ends_by_redundant_and_support, subset_check_exons

JUNCS = (Junc(200, 300), Junc(400, 500))


def _variant(start, end, nreads, extra_read=None):
    """an end variant with representative ends start and end, its reads there, and
    optionally one more read with other ends"""
    iso = Isoform('chr1', '+', JUNCS, start, end)
    iso.reads.extend([SimpleNamespace(start=start, end=end) for _ in range(nreads)])
    if extra_read is not None:
        iso.reads.append(SimpleNamespace(start=extra_read[0], end=extra_read[1]))
    return iso


def test_normalize_ends_one_isoform_with_furthest_ends():
    # the furthest start and the furthest end come from different variants
    variants = [_variant(100, 600, 5), _variant(50, 550, 2), _variant(120, 700, 1)]
    for max_ends in (1, 3):
        kept = filter_ends_by_redundant_and_support([_variant(v.start, v.end, len(v.reads)) for v in variants],
                                                    sjc_support=1, se_support=3, max_ends=max_ends, normalize_ends=True)
        assert len(kept) == 1
        assert (kept[0].start, kept[0].end, kept[0].num_reads) == (50, 700, 8)


def test_normalize_ends_furthest_read_not_representative():
    # a read past its variant's representative ends sets the start
    variants = [_variant(100, 600, 5), _variant(50, 550, 1, extra_read=(30, 540))]
    kept = filter_ends_by_redundant_and_support(variants, sjc_support=1, se_support=3, max_ends=1,
                                                normalize_ends=True)
    assert (kept[0].start, kept[0].end, kept[0].num_reads) == (30, 600, 7)
    # the subset check still uses the best variant's representative ends
    exons = subset_check_exons(kept[0])
    assert (exons[0].start, exons[-1].end) == (100, 600)


def test_normalize_ends_support_is_the_chain_total():
    kept = filter_ends_by_redundant_and_support([_variant(100, 600, 1), _variant(50, 550, 1)],
                                                sjc_support=3, se_support=3, max_ends=1, normalize_ends=True)
    assert kept == []


def test_without_normalize_ends_max_ends_still_applies():
    kept = filter_ends_by_redundant_and_support([_variant(100, 600, 5), _variant(50, 550, 4)],
                                                sjc_support=2, se_support=3, max_ends=2, normalize_ends=False)
    assert len(kept) == 2


def test_subset_check_uses_the_ends_chosen_without_normalize_ends():
    variants = [(100, 600, 5), (50, 550, 2), (120, 700, 1)]
    results = {}
    for normalize_ends in (False, True):
        kept = filter_ends_by_redundant_and_support([_variant(*v) for v in variants], sjc_support=1, se_support=3,
                                                    max_ends=1, normalize_ends=normalize_ends)
        results[normalize_ends] = subset_check_exons(kept[0])
    assert results[True] == results[False]
    assert (results[True][0].start, results[True][-1].end) == (100, 600)
