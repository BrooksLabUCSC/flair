from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import (BEST_END_WINDOW, densest_end, normalize_chain_ends,
                                       filter_ends_by_redundant_and_support, subset_check_exons)

JUNCS = (Junc(200, 300), Junc(400, 500))


def _chain(read_ends):
    iso = Isoform('chr1', '+', JUNCS)
    iso.reads.extend(SimpleNamespace(start=s, end=e) for s, e in read_ends)
    return iso


def test_densest_end_prefers_a_cluster_over_scattered_truncated_reads():
    # 3 reads at the TSS, 5 truncated reads spread out further in
    starts = [100, 102, 105] + [140, 170, 200, 230, 260]
    assert densest_end(starts, BEST_END_WINDOW, min) == 100


def test_densest_end_ties_go_to_the_outer_end():
    assert densest_end([100, 200], BEST_END_WINDOW, min) == 100
    assert densest_end([600, 700], BEST_END_WINDOW, max) == 700
    assert densest_end([150], BEST_END_WINDOW, min) == 150


def test_normalize_chain_ends_furthest_and_best_supported():
    iso = normalize_chain_ends(_chain([(100, 600), (102, 605), (104, 598), (30, 590), (180, 700)]))
    assert (iso.start, iso.end) == (30, 700)
    assert iso.best_ends == (100, 605)
    exons = subset_check_exons(iso)
    assert (exons[0].start, exons[-1].end) == (100, 605)


def test_without_normalize_ends_max_ends_still_applies():
    variants = []
    for start, end, nreads in ((100, 600, 5), (50, 550, 4)):
        iso = Isoform('chr1', '+', JUNCS, start, end)
        iso.reads.extend(SimpleNamespace(start=start, end=end) for _ in range(nreads))
        variants.append(iso)
    kept = filter_ends_by_redundant_and_support(variants, sjc_support=2, se_support=3, max_ends=2)
    assert len(kept) == 2
