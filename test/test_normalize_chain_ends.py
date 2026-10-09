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
    assert densest_end(starts, BEST_END_WINDOW, min) == 102


def test_densest_end_in_a_cluster_is_its_mode_not_its_outer_edge():
    # every position in the cluster has all its reads within the window
    assert densest_end([990, 1000, 1000, 1000, 1012, 1020], BEST_END_WINDOW, max) == 1000
    assert densest_end([990, 1000, 1012, 1020], BEST_END_WINDOW, max) == 1012


def test_densest_end_ties_go_to_the_middle_then_the_outer_end():
    assert densest_end([100, 110, 120], BEST_END_WINDOW, min) == 110
    assert densest_end([100, 110, 120, 130], BEST_END_WINDOW, max) == 120


def test_densest_end_of_scattered_read_ends_is_the_outer_one():
    # no read end has another within the window: there is no cluster
    assert densest_end([100, 300, 500], BEST_END_WINDOW, min) == 100
    assert densest_end([100, 300, 500], BEST_END_WINDOW, max) == 500
    assert densest_end([100, 200], BEST_END_WINDOW, min) == 100
    assert densest_end([600, 700], BEST_END_WINDOW, max) == 700
    assert densest_end([150], BEST_END_WINDOW, min) == 150


def test_normalize_chain_ends_furthest_and_best_supported():
    iso = normalize_chain_ends(_chain([(100, 600), (102, 605), (104, 598), (30, 590), (180, 700)]))
    assert (iso.start, iso.end) == (30, 700)
    assert iso.best_ends == (102, 600) and iso.three_prime_clustered
    exons = subset_check_exons(iso)
    assert (exons[0].start, exons[-1].end) == (102, 600)


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


def test_normalize_chain_ends_records_a_scattered_3_prime_end():
    # on +, the 3' end is the read ends; they are scattered, the starts not
    iso = normalize_chain_ends(_chain([(100, 600), (102, 900), (104, 1300)]))
    assert iso.best_ends == (102, 1300) and not iso.three_prime_clustered
