from types import SimpleNamespace
from flair.isoform_data import Exon, Isoform, Junc
from flair.flair_transcriptome import (CandidateIsoforms, _longest_supported_read_ends, _single_exon_end_variants,
                                       filter_firstpass_isos, group_se_by_overlap, normalize_chain_ends)


def _read(start, end, strand, polyA=(0, 0), intprim=(False, False)):
    return SimpleNamespace(start=start, end=end, strand=strand, polyA=polyA, intprim=intprim)


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


def _spliced(name, strand, exons, read_ends=None):
    "a spliced candidate, its ends set from its reads' as the first pass sets them"
    reads = [SimpleNamespace(start=s, end=e) for s, e in (read_ends or [(exons[0][0], exons[-1][1])] * 5)]
    iso = Isoform('chr1', strand, tuple(Junc(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)), reads=reads)
    iso.name = name
    return normalize_chain_ends(iso)


def _annots(strand, exons, single=()):
    "an annotation with one spliced transcript, and single-exon ones"
    return SimpleNamespace(transcripts=[('tx', 'gene', strand)],
                           transcript_to_exons={('tx', 'gene'): [Exon(s, e) for s, e in exons]},
                           all_annot_SE={'+': [], '-': [], strand: [Exon(s, e, 'single_gene') for s, e in single]})


COMPREHENSIVE = SimpleNamespace(filter='comprehensive', sjc_support=1, end_window=100, ss_window=10)


def test_single_exon_isoform_starting_in_a_3_prime_terminal_exon_is_a_utr_fragment():
    # starting in the spliced isoform's last exon and running past its end
    spliced = _spliced('spliced', '+', [(1000, 2000), (3000, 4000)])
    fragment, intronic = _se('fragment', 3500, 6000, '+', 5), _se('intronic', 2200, 2800, '+', 5)
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([spliced, fragment, intronic]), None, {})
    assert sorted(firstpass) == ['intronic', 'spliced']


def test_the_5_prime_terminal_exon_does_not_make_a_fragment():
    spliced = _spliced('spliced', '+', [(1000, 2000), (3000, 4000)])
    upstream = _se('upstream', 1500, 2500, '+', 5)  # starts in the first exon, runs into the intron
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([spliced, upstream]), None, {})
    assert sorted(firstpass) == ['spliced', 'upstream']


def test_an_annotated_3_prime_terminal_exon_makes_a_fragment_on_minus():
    # on -, the 3' terminal exon is the first in genomic order, and the 5' end the right end
    annots = _annots('-', [(1000, 5000), (6000, 7000)])
    fragment, other = _se('fragment', 1500, 4000, '-', 5), _se('other', 500, 900, '-', 5)
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([fragment, other]), annots, {})
    assert sorted(firstpass) == ['other']


def test_a_fragment_starting_as_far_past_the_terminal_exon_as_a_read_end_cluster_reaches():
    # 5' ends 13 and 40 bp past the 3' terminal exon's end at 5000; BEST_END_WINDOW is 25
    annots = _annots('-', [(1000, 5000), (6000, 7000)])
    near, past = _se('near', 1500, 5013, '-', 5), _se('past', 1500, 5040, '-', 5)
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([near]), annots, {})
    assert sorted(firstpass) == []
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([past]), annots, {})
    assert sorted(firstpass) == ['past']


def test_a_single_exon_isoform_matching_an_annotated_single_exon_transcript_is_not_a_fragment():
    annots = _annots('-', [(1000, 5000), (6000, 7000)], single=[(1500, 4000)])
    annotated, fragment = _se('annotated', 1510, 3990, '-', 5), _se('fragment', 3000, 4500, '-', 5)
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([annotated, fragment]), annots, {})
    assert sorted(firstpass) == ['annotated']


def test_a_spliced_isoform_whose_3_prime_end_is_not_a_cluster_makes_no_fragment():
    # two reads of the chain, their 3' ends far apart: the outer one, 5000, is its end
    scattered = _spliced('scattered', '+', [(1000, 2000), (3000, 5000)], read_ends=[(1000, 5000), (1000, 3500)])
    single = _se('single', 4000, 6000, '+', 5)
    firstpass, _ = filter_firstpass_isos(COMPREHENSIVE, _candidates([scattered, single]), None, {})
    assert sorted(firstpass) == ['scattered', 'single']


def test_internally_primed_reads_are_left_out_before_stranding():
    # two reads internally primed at their right ends
    reads = [_read(1000, 2000, '+', (0, 20)), _read(1010, 1990, '+', (0, 15)), _read(1005, 1995, '+', (0, 18)),
             _read(1002, 1998, '+', intprim=(False, True)), _read(1003, 1997, '+', intprim=(False, True))]
    ((strand, reads_kept),) = _groups(reads, trust_strand=False, se_support=3)
    assert strand == '+' and reads_kept == [(1000, 2000), (1005, 1995), (1010, 1990)]


def test_with_trust_strand_only_the_3_prime_end_is_checked():
    # priming at the start of a + read is at its 5' end: not priming
    five_prime = [_read(1000, 2000, '+', intprim=(True, False))] * 2
    three_prime = [_read(1000, 2000, '+', intprim=(False, True))] * 2
    assert _groups(five_prime, trust_strand=True) == [('+', [(1000, 2000), (1000, 2000)])]
    assert _groups(three_prime, trust_strand=True) == []


def test_total_rna_keeps_internally_primed_reads():
    reads = [_read(1000, 2000, '+', intprim=(False, True))] * 3
    isoform = Isoform('chr1', '+', (), reads=reads)
    assert list(group_se_by_overlap('chr1', isoform, 3, True)) == []
    ((_, strand, kept),) = group_se_by_overlap('chr1', isoform, 3, True, total_rna=True)
    assert strand == '+' and len(kept) == 3


def _variants(read_ends, support=3, min_frac=0.05):
    reads = [SimpleNamespace(start=s, end=e) for s, e in read_ends]
    variants = _single_exon_end_variants(Isoform('chr1', '+', ()), reads, support=support, min_frac=min_frac)
    return sorted((iso.start, iso.end, len(iso_reads)) for iso, iso_reads in variants)


def test_truncated_reads_of_a_long_single_exon_transcript_are_one_isoform():
    # 4 full-length reads, and 10 truncated at scattered starts: no other start cluster
    reads = [(1000, 10000)] * 4 + [(2000 + 700 * i, 10000) for i in range(10)]
    assert _variants(reads) == [(1000, 10000, 14)]


def test_single_exon_isoform_at_each_clustered_pair_of_ends():
    reads = [(1000, 2000)] * 5 + [(1000, 2600)] * 5 + [(1010, 2300)]
    assert _variants(reads) == [(1000, 2000, 6), (1000, 2600, 5)]


def test_without_clustered_ends_one_isoform_at_the_longest_supported_read():
    reads = [(1000, 2000), (1010, 2010), (1500, 3000), (1800, 2500)]
    assert _variants(reads) == [(1000, 2000, 4)]


def test_a_chance_cluster_at_a_deeply_covered_locus_is_no_end():
    # 3 reads are a cluster by count, not 5% of 203
    reads = [(1000, 5000)] * 200 + [(3000, 5000)] * 3
    assert _variants(reads) == [(1000, 5000, 203)]


def test_longest_supported_ends_beside_fragment_clusters():
    # 5' and 3' fragments cluster; the 3 full-length reads are too few for 5% of 63
    reads = [(1000, 1500)] * 30 + [(2500, 3000)] * 30 + [(1000, 3000)] * 2 + [(1005, 2995)]
    assert _variants(reads) == [(1000, 1500, 30), (1000, 3000, 3), (2500, 3000, 30)]
    variants = _single_exon_end_variants(Isoform('chr1', '+', ()), [SimpleNamespace(start=s, end=e) for s, e in reads],
                                         support=3, min_frac=0.05)
    assert [(iso.start, iso.longest_supported) for iso, _ in variants] == [(1000, False), (1000, True), (2500, False)]


def _longest_flags(read_ends, support=3):
    reads = [SimpleNamespace(start=s, end=e) for s, e in read_ends]
    variants = _single_exon_end_variants(Isoform('chr1', '+', ()), reads, support=support, min_frac=0.05)
    return [(iso.start, iso.end, iso.longest_supported) for iso, _ in variants]


def test_longest_supported_ends_outside_any_cluster_defer_to_the_clustered_pair():
    # 2 reads overhang the cluster's 5' end by 100 bp: an outlier, not a start
    reads = [(1000, 2000)] * 5 + [(900, 2010)] * 2
    assert _variants(reads) == [(1000, 2000, 7)]
    assert _longest_flags(reads) == [(1000, 2000, True)]


def test_clustered_pair_takes_the_longest_supported_ends_in_a_cluster():
    # the end at 2000 is a cluster of 5, with the 3' fragments, but only 3 reads
    # have both ends in (1000, 2000), fewer than support 4
    reads = [(1000, 1900)] * 10 + [(1000, 2000)] * 3 + [(1300, 2000)] * 2
    assert _variants(reads, support=4) == [(1000, 2000, 15)]
    assert _longest_flags(reads, support=4) == [(1000, 2000, True)]
