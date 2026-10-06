from types import SimpleNamespace
from flair.isoform_data import Isoform, Junc
from flair.flair_transcriptome import (CandidateIsoforms, _longest_supported_read_ends, filter_firstpass_isos,
                                       group_se_by_overlap)


def _read(start, end, strand, polyA=(0, 0)):
    return SimpleNamespace(start=start, end=end, strand=strand, polyA=polyA)


def _groups(reads, trust_strand, se_support=2):
    isoform = Isoform('chr1', '+', (), reads=reads)
    return sorted((strand, [(r.start, r.end) for r in group])
                  for _, strand, group in group_se_by_overlap('chr1', isoform, se_support, trust_strand))


def test_trust_strand_groups_each_strand_separately():
    # a sense and an antisense transcript at the same locus
    reads = [_read(1000, 2000, '+'), _read(1010, 1990, '+'), _read(1005, 1995, '-')]
    assert _groups(reads, trust_strand=True) == [('+', [(1000, 2000), (1010, 1990)]), ('-', [(1005, 1995)])]


def test_without_trust_strand_overlapping_reads_group_and_poly_a_gives_the_strand():
    # alignment strands mean nothing; two right poly(A) tails make the group +
    reads = [_read(1000, 2000, '-', (0, 20)), _read(1010, 1990, '+', (0, 15)), _read(1005, 1995, '-')]
    assert _groups(reads, trust_strand=False) == [('+', [(1000, 2000), (1005, 1995), (1010, 1990)])]


def test_without_a_poly_a_consensus_the_group_is_dropped():
    reads = [_read(1000, 2000, '+', (0, 20)), _read(1010, 1990, '+', (20, 0))]
    assert _groups(reads, trust_strand=False) == []


def _candidates(isoforms):
    candidates = CandidateIsoforms()
    for iso in isoforms:
        candidates.add(iso)
    return candidates


def _se(name, start, end, strand, nreads, cluster=None):
    iso = Isoform('chr1', strand, (), start, end, reads=[SimpleNamespace(start=start, end=end)] * nreads)
    iso.name = name
    iso.end_variant_cluster = cluster
    return iso


def test_single_exon_filter_only_compares_isoforms_on_its_strand():
    spliced = Isoform('chr1', '+', (Junc(2000, 3000),), 1000, 4000, reads=[SimpleNamespace(start=1000, end=4000)] * 5)
    spliced.name = 'spliced'
    # each inside the spliced isoform's first exon, and inside the other's span
    sense, antisense = _se('sense', 1100, 1900, '+', 5), _se('antisense', 1100, 1900, '-', 5)
    args = SimpleNamespace(filter='comprehensive', sjc_support=1, end_window=100, ss_window=10)
    firstpass, _ = filter_firstpass_isos(args, _candidates([spliced, sense, antisense]), None, {})
    assert sorted(firstpass) == ['antisense', 'spliced']


def test_end_variants_of_one_cluster_are_not_compared():
    args = SimpleNamespace(filter='nosubset', sjc_support=1, end_window=100, ss_window=10)
    same = [_se('long', 1000, 3000, '+', 5, cluster='c'), _se('short', 1100, 1900, '+', 5, cluster='c')]
    firstpass, _ = filter_firstpass_isos(args, _candidates(same), None, {})
    assert sorted(firstpass) == ['long', 'short']
    other = [_se('long', 1000, 3000, '+', 5, cluster='c'), _se('short', 1100, 1900, '+', 5, cluster='d')]
    firstpass, _ = filter_firstpass_isos(args, _candidates(other), None, {})
    assert sorted(firstpass) == ['long']


SEG1_READS = [(32186476, 32189820), (32186476, 32189823), (32186479, 32188247), (32186479, 32188247), (32186480, 32188145),
              (32186485, 32188140), (32186486, 32188134), (32186515, 32188421), (32187393, 32188140), (32187395, 32189017),
              (32187395, 32189877), (32187794, 32190167), (32188011, 32188426), (32188011, 32189054), (32188011, 32189241),
              (32188013, 32188303), (32188013, 32188336), (32188013, 32188476), (32188013, 32188700), (32188013, 32188752),
              (32188013, 32188876), (32188020, 32188303), (32188022, 32188814)]


def _ends(read_ends):
    return _longest_supported_read_ends([SimpleNamespace(start=s, end=e) for s, e in read_ends])


def test_single_exon_ends_are_the_longest_read_whose_both_ends_have_support():
    # the reads running furthest, to 32189877 and 32190167, have no other read ending near them
    assert _ends(SEG1_READS) == (32186476, 32189823)


def test_a_read_sticking_out_past_the_others_does_not_set_the_ends():
    assert _ends([(1000, 2000), (1010, 1990), (500, 2000), (1005, 2600)]) == (1000, 2000)


def test_without_a_read_supported_at_both_ends_the_longest_read_is_used():
    assert _ends([(1000, 2000), (1500, 2600), (1800, 3000)]) == (1800, 3000)


def test_the_fallback_skips_reads_more_than_twice_as_long_as_the_next():
    # 5000 is more than twice 2000, which isn't more than twice 1100
    assert _ends([(0, 5000), (1000, 3000), (1300, 2400), (1500, 2100)]) == (1000, 3000)
