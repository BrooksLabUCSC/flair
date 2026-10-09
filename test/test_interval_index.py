import random
from flair.interval_index import IntervalIndex
from flair.intron_support import IntronSupport


def _brute_overlap(intervals, start, end, slack=0):
    return [d for s, e, d in intervals if max(s, start - slack) < min(e, end + slack)]


def test_overlap_matches_scan_in_insertion_order():
    rng = random.Random(1)
    intervals = []
    index = IntervalIndex()
    for i in range(500):
        s = rng.randrange(0, 100000)
        e = s + rng.choice((0, 1, rng.randrange(1, 300), rng.randrange(1, 20000)))
        intervals.append((s, e, i))
        index.add(s, e, i)
    for _ in range(300):
        s = rng.randrange(-1000, 101000)
        e = s + rng.randrange(0, 5000)
        slack = rng.choice((0, 0, 50))
        assert index.overlap(s, e, slack) == _brute_overlap(intervals, s, e, slack)


def test_add_after_query_rebuilds():
    index = IntervalIndex()
    index.add(100, 200, 'a')
    assert index.overlap(150, 160) == ['a']
    index.add(10, 1000, 'b')
    assert index.overlap(150, 160) == ['a', 'b']
    assert index.overlap(200, 300) == ['b']
    assert IntervalIndex().overlap(0, 10) == []


def test_intron_subset_matches_scan():
    rng = random.Random(2)
    support = IntronSupport()
    for _ in range(400):
        start = rng.randrange(0, 200000)
        end = start + rng.randrange(support.min_intron_size, 50000)
        support.add_support('chr1', start, end, rng.choice('+-'), rng.randrange(1, 5))
    support.build_indexes()
    for _ in range(100):
        start = rng.randrange(0, 250000)
        end = start + rng.randrange(1, 30000)
        expect = sorted((i.start, i.end, i.strand, i.read_support_cnt) for i in support.introns('chr1')
                        if i.end > start and i.start < end)
        sub = support.subset_for_region('chr1', start, end)
        assert sorted((i.start, i.end, i.strand, i.read_support_cnt) for i in sub.introns('chr1')) == expect
