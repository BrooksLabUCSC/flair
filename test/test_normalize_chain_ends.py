from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import (BEST_END_WINDOW, _end_variants, densest_end, normalize_chain_ends,
                                       subset_check_exons)

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


def _reads(*groups):
    return [SimpleNamespace(start=start, end=end) for start, end, nreads in groups for _ in range(nreads)]


def test_end_groups_with_support_are_kept_up_to_max_ends():
    reads = _reads((100, 600, 5), (0, 450, 4), (150, 850, 1))
    variants = _end_variants(_chain([]), reads, max_ends=2, end_window=20, support=2)
    assert sorted((v.start, v.end, len(variant_reads)) for v, variant_reads in variants) == [(0, 450, 4), (100, 600, 6)]


def test_fewer_than_two_supported_end_groups_give_one_isoform_at_the_densest_ends():
    reads = _reads((100, 600, 5), (0, 450, 1))
    ((iso, iso_reads),) = _end_variants(_chain([]), reads, max_ends=2, end_window=20, support=2)
    assert (iso.start, iso.end, len(iso_reads)) == (100, 600, 6)
