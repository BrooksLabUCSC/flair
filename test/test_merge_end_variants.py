from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import filter_ends_by_redundant_and_support

JUNCS = (Junc(200, 300), Junc(400, 500))


def _variant(start, end, nreads):
    iso = Isoform('chr1', '+', JUNCS, start, end)
    iso.reads.extend([object()] * nreads)
    return iso


def test_normalize_ends_one_isoform_with_furthest_ends():
    # the furthest start and the furthest end come from different variants
    variants = [_variant(100, 600, 5), _variant(50, 550, 2), _variant(120, 700, 1)]
    for max_ends in (1, 3):
        kept = filter_ends_by_redundant_and_support([_variant(v.start, v.end, len(v.reads)) for v in variants],
                                                    sjc_support=1, se_support=3, max_ends=max_ends, normalize_ends=True)
        assert len(kept) == 1
        assert (kept[0].start, kept[0].end, kept[0].num_reads) == (50, 700, 8)


def test_normalize_ends_support_is_the_chain_total():
    kept = filter_ends_by_redundant_and_support([_variant(100, 600, 1), _variant(50, 550, 1)],
                                                sjc_support=3, se_support=3, max_ends=1, normalize_ends=True)
    assert kept == []


def test_without_normalize_ends_max_ends_still_applies():
    kept = filter_ends_by_redundant_and_support([_variant(100, 600, 5), _variant(50, 550, 4)],
                                                sjc_support=2, se_support=3, max_ends=2, normalize_ends=False)
    assert len(kept) == 2
